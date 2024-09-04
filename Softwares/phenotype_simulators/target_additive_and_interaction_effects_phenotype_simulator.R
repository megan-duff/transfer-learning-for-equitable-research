#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}

# Check if genio is installed; if not, install it
if (!requireNamespace("genio", quietly = TRUE)) {
  install.packages("genio", repos = "https://cloud.r-project.org/")
}

library(optparse)
library(genio) 

# Define command-line arguments
option_list <- list(
  make_option(c("-p", "--plink_file"), type="character", default=NULL, help="Path to PLINK files", metavar="FILE"),
  make_option(c("-a", "--af_file"), type="character", default=NULL, help="Path to allele frequency file", metavar="FILE"),
  make_option(c("-c", "--source_causal_snp_file"), type="character", default=NULL, help="Path to causal loci file", metavar="FILE"),
  make_option(c("-l", "--h2"), type="numeric", default=NULL, help="Overall Heritability", metavar="NUM"),
  make_option(c("-e", "--e"), type="numeric", default=NULL, help="Percentage of heritability attributed to interaction effects", metavar="NUM"),
  make_option(c("-n", "--num_interaction_terms"), type="integer", default=NULL, help="Number of interaction terms", metavar="NUM"),
  make_option("--source_par_file", type="character", default=NULL, help="Path to source par file", metavar="FILE"),
  make_option("--source_interaction_file", type="character", default=NULL, help="Path to source interaction file", metavar="FILE"),
  make_option(c("-g", "--genetic_correlation"), type="integer", default=NULL, help="correlation between source and target set effects", metavar="NUM"),
  make_option(c("-o", "--output"), type="character", default="output", help="Output file prefix", metavar="PREFIX")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Extract arguments
plink_file <- args$plink_file
af_file <- args$af_file
causal_snp_file <- args$causal_snp_file
h2 <- args$h2
e <- args$e
num_interaction_terms <- args$num_interaction_terms
output_file_prefix <- args$output

source_causal_loci_path<- args$source_causal_snp_file 
source_par_file_path<- args$source_par_file 
source_interaction_file_path<- args$source_interaction_file
genetic_correlation<- args$genetic_correlation


# Load in data function
load_in_data <- function(file) {
  plink_data <- read_plink(file)
  print("Successfully read in .bed/.bim/.fam files!")
  return(plink_data)
}

# Get bim function
get_bim <- function(plink_data) {
  bim <- plink_data$bim
  return(bim)
}

# SD data function
sd_data <- function(plink_data, af_file) {
  af <- read.table(af_file, header = FALSE)
  print("Read in AF file")
  colnames(af) <- c("#CHROM", "ID", "REF", "ALT", "PROVISIONAL_REF?", "ALT_FREQS", "OBS_CT")
  genotypes <- plink_data$X
  print("Grabbed genotypes")
  sd_data <- (genotypes - 2*af$ALT_FREQS)/sqrt(2*af$ALT_FREQS*(1-af$ALT_FREQS))
  print("Standardized!")
  return(sd_data)
}

target_non_linear_pheno_sim <- function(plink_data, sd_data, source_causal_loci, source_par_file, source_interaction_file, genetic_correlation, h2, e, num_interaction_terms, output_file_prefix) {
  
  # Get genotype matrix - rows are SNPs and columns as samples 
  genotype <- sd_data
  sample_size <- dim(genotype)[2]
  
  # Read in causal loci file 
  causal_loci_vector <- read.table(source_causal_loci)
  num_of_causal_snps <- dim(causal_loci_vector)[1]
  
  print(paste("Number of causal SNPs:", num_of_causal_snps, sep=" "))

  # Add quotes around each rs number for the causal SNP list to match format 
  causal_genotype_quote <- sapply(causal_loci_vector, function(x) sprintf("%s", x))
  
  # Subset genoype matrix to only include causal loci 
  genotype_matrix_causal <- genotype[rownames(genotype) %in% causal_genotype_quote, ]
  
  print(paste("Dimension of causal genotype matrix:", dim(genotype_matrix_causal), sep=" "))
  
  # Get reference allele from bim file of plink_data 
  bim <- plink_data$bim
  causal_bim <- bim[bim$id %in% causal_genotype_quote, ]
  ref_alleles <- causal_bim$ref
  
  interaction_snps <- matrix(nrow=0, ncol=2)
  interaction_matrix <- matrix(nrow=0, ncol=dim(genotype_matrix_causal)[2])
  
  # Read in interaction snps file to create target interaction matrix 
  interaction_file <- read.table(source_interaction_file, sep = "\t", header=FALSE)
  
  # Select interaction terms and create vector of interaction terms for each individual
  for (i in 1:num_interaction_terms){
    interaction_term <- c(match(interaction_file[i,1], rownames(genotype_matrix_causal)),
                          match(interaction_file[i,2], rownames(genotype_matrix_causal)))
    
    interaction_vector <- genotype_matrix_causal[interaction_term[1],] * genotype_matrix_causal[interaction_term[2],]
    interaction_matrix <- rbind(interaction_matrix, interaction_vector)
    interaction_snps <- rbind(interaction_snps, c(rownames(genotype_matrix_causal)[interaction_term[1]], rownames(genotype_matrix_causal)[interaction_term[2]]))
  }
  
  print("Starting beta simulation for causal SNPs and interaction terms...")
  
  h2_interaction <- h2*e
  h2_additive <- h2 - h2_interaction
  
  # Read in source additive betas file and extract additive beta values
  par_file <- read.table(source_par_file, sep = "\t", header=TRUE)
  source_additive_betas <- as.numeric(par_file$Scaled.Additive.Effects)
  
  # Read in source interactions betas file and extract interactions beta values
  interaction_file <- read.table(source_interaction_file, sep = "\t", header=FALSE)
  source_interaction_betas <- as.numeric(interaction_file[,3])
  
  # Draw beta values for additive SNPs
  additive_betas <- rnorm(num_of_causal_snps, mean = source_additive_betas, sd = sqrt((1-genetic_correlation)*h2_additive/num_of_causal_snps))
  interaction_betas <- rnorm(num_interaction_terms, mean = interaction_file[,3], sd = sqrt((1-genetic_correlation)*h2_interaction/num_interaction_terms))
  
  # Compute allele frequencies 
  allele_frequencies <- rowSums(genotype_matrix_causal) / (2 * ncol(genotype_matrix_causal))
  
  print("Starting phenotype simulation for each individual...")

  # Draw error values 
  error_vector <- rnorm(sample_size, mean = 0, sd = sqrt(1-h2))
  
  # Compute XB for each individual
  additive_beta_matrix <- matrix(additive_betas, nrow = length(additive_betas), ncol = 1)
  additive_genetic_contribution_vector <- t(genotype_matrix_causal) %*% additive_beta_matrix
  
  # Compute ZY for each individual
  interaction_beta_matrix <- matrix(interaction_betas, nrow = length(interaction_betas), ncol = 1)
  interaction_genetic_contribution_vector <- t(interaction_matrix) %*% interaction_beta_matrix

  # Compute XB + ZY
  genetic_contribution_vector <- interaction_genetic_contribution_vector + additive_genetic_contribution_vector
  
  # Compute phenotype
  phenotype_vector <- genetic_contribution_vector + error_vector
  
  # Compute Var(XB), Var(ZY), Var(e)
  var_XB = var(additive_genetic_contribution_vector)
  var_ZY = var(interaction_genetic_contribution_vector)
  var_error = var(error_vector)
  var = var(genetic_contribution_vector)
  
  # Print relevant statistics for phenotype simulation 
  print(paste("Additive Heritability:", var_XB/var(phenotype_vector), sep=" "))
  print(paste("Interaction Heritability:", var_ZY/var(phenotype_vector), sep=" "))
  print(paste("Error % of Variance:", var(error_vector)/var(phenotype_vector), sep=" "))
  print(paste("Heritability:", var/var(phenotype_vector), sep=" "))
  
  # Create error file 
  error_file <- cbind(plink_data$fam$fam, plink_data$fam$id, additive_genetic_contribution_vector, interaction_genetic_contribution_vector, genetic_contribution_vector, error_vector)
  colnames(error_file) <- c("id_1","id_2", "XB", "ZY", "Genetic_Contribution", "Error")
  
  # Create phenotype file 
  phenotype_file <- cbind(plink_data$fam$fam, plink_data$fam$id, phenotype_vector)
  
  # Create interaction file 
  interaction_file <- cbind(interaction_snps, interaction_beta_matrix)
  
  # Create par_file 
  par_file <- cbind(rownames(genotype_matrix_causal), ref_alleles, allele_frequencies, additive_betas)
  colnames(par_file) <- c("QTL", "RefAllele", "Frequency", "Additive_Effects")

  # Write phenotype file
  print("Writing phenotype file...")
  write.table(phenotype_file, file = paste(output_file_prefix, ".phen", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=FALSE) 
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".phen", sep = ""), sep=" "))
  
  # Write interaction file
  print("Writing interaction file...")
  write.table(interaction_file, file = paste(output_file_prefix, "_interactions.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE) 
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, "_interactions.txt", sep = ""), sep=" "))
  
  # Write par file
  print("Writing par file...")
  write.table(par_file, file = paste(output_file_prefix, ".par", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=TRUE)   
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".par", sep = ""), sep=" "))

  # Write error file 
  print("Writing error file...")
  write.table(error_file, file = paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=TRUE)
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep=" "))
  
}

plink_data<-load_in_data(plink_file) 

sd_data<-sd_data(plink_data, af_file) 

target_non_linear_pheno_sim(
                        sd_data=sd_data,
                        plink_data=plink_data,
                        source_causal_loci=source_causal_loci_path, 
                        source_par_file=source_par_file_path, 
                        source_interaction_file=source_interaction_file_path,
                        genetic_correlation=genetic_correlation, 
                        h2=h2, 
                        e=e,
                        num_interaction_terms=num_interaction_terms,
                        output_file_prefix=output_file_prefix)
