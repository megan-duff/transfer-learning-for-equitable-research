#!/usr/bin/env Rscript

load_in_data <- function(file) {
  # Read in PLINK fileset
  plink_data <- read_plink(file)
  print("Successfully read in .bed/.bim/.fam files!")
  return(plink_data)
}

get_bim <- function(plink_data) {
  bim <- plink_data$bim
  return(bim)
}

sd_data <- function(plink_data, af_file) {
  af <- read.table(af_file, header = FALSE)
  print("read in AF file")
  colnames(af)<-c("#CHROM",	"ID",	"REF",	"ALT",	"ALT_FREQS",	"OBS_CT")
  genotypes <- plink_data$X
  print("Grabbed genotypes")
  sd_data <- (genotypes - 2*af$ALT_FREQS)/sqrt(2*af$ALT_FREQS*(1-af$ALT_FREQS))
  print("standardized!")
  return(sd_data)
}

source_linear_pheno_sim <- function(plink_data, sd_data, causal_loci, h2, output_file_prefix) {
  
  # Get genotype matrix - rows are SNPs and columns as samples 
  genotype <- sd_data
  sample_size <- dim(genotype)[2]
  
  # Read in causal loci file 
  causal_loci_vector <- read.table(causal_loci)
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
  
  print("Starting beta simulation for causal SNPs...")
  
  # Draw beta values 
  betas <- rnorm(num_of_causal_snps, mean = 0, sd = sqrt(h2/num_of_causal_snps))
  
  print(paste("Mean of betas:", mean(betas), sep=" "))
  
  # Compute allele frequencies 
  allele_frequencies <- rowSums(genotype_matrix_causal) / (2 * ncol(genotype_matrix_causal))
  
  # Create par_file with relevant information 
  par_file <- cbind(rownames(genotype_matrix_causal), ref_alleles, allele_frequencies, betas)
  colnames(par_file) <- c("QTL", "RefAllele", "Frequency", "Effect")
  
  print("Starting phenotype simulation for each individual...")
  
  # Compute XB for each individual
  beta_matrix <- matrix(betas, nrow = length(betas), ncol = 1)
  
  print(paste("Dimension of t(genotype_matrix_causal):", dim(t(genotype_matrix_causal))[1],dim(t(genotype_matrix_causal))[2], sep=" "))
  
  print(paste("Dimension of beta matrix:", dim(beta_matrix)[1],dim(beta_matrix)[2], sep=" "))
  
  genetic_contribution_vector <- t(genotype_matrix_causal) %*% beta_matrix
  
  print(paste("Dimension of genetic_contribution_vector:", dim(genetic_contribution_vector)[1], dim(genetic_contribution_vector)[2], sep=" "))
    
  # Compute Var(XB)
  var = var(genetic_contribution_vector)
  
  # Draw error values 
  error_vector <- rnorm(sample_size, mean = 0, sd = sqrt(1-h2))
  
  # Compute phenotype
  phenotype_vector <- genetic_contribution_vector + error_vector
  
  print(paste("Range of phenotype vector:", range(phenotype_vector)[1], range(phenotype_vector)[2], sep = " "))
  
  # Print relevant statistics for phenotype simulation 
  print(paste("Heritability:", var/var(phenotype_vector), sep=" "))
  
  # Create error file 
  error_file <- cbind(plink_data$fam$fam, plink_data$fam$id, genetic_contribution_vector, error_vector)
  colnames(error_file) <- c("id_1","id_2", "XB", "Error")
  
  # Create phenotype file 
  phenotype_file <- cbind(plink_data$fam$fam, plink_data$fam$id, phenotype_vector)

  # Write phenotype file
  print("Writing phenotype file...")
  write.table(phenotype_file, file = paste(output_file_prefix, ".phen", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=FALSE) 
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".phen", sep = ""), sep=" "))
  
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

source_linear_pheno_sim_subset <- function(plink_data, sd_data, source_causal_loci, source_par_file, h2,  output_file_prefix) {

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
  
  print(paste("Number of causal SNPs:", dim(genotype_matrix_causal)[1], sep=" "))
  print(paste("Sample Size:", dim(genotype_matrix_causal)[2], sep=" "))
  
  # Get reference allele from bim file of plink_data 
  bim <- plink_data$bim
  causal_bim <- bim[bim$id %in% causal_genotype_quote, ]
  ref_alleles <- causal_bim$ref
  
  # Read in source betas file and extract beta values
  par_file <- read.table(source_par_file, sep = "\t", header=TRUE)
  source_betas <- as.numeric(par_file$Effect)
  
  print("Starting beta simulation for causal SNPs...")
  
  # Draw beta values 
  betas <- source_betas
  
  # Compute allele frequencies 
  allele_frequencies <- rowSums(genotype_matrix_causal) / (2 * ncol(genotype_matrix_causal))
  
  # Create par_file with relevant information 
  par_file <- cbind(rownames(genotype_matrix_causal), ref_alleles, allele_frequencies, betas)
  colnames(par_file) <- c("QTL", "RefAllele", "Frequency", "Effect")
  
  print("Starting phenotype simulation for each individual...")
  
  # Compute XB for each individual
  beta_matrix <- matrix(betas, nrow = length(betas), ncol = 1)
  
  genetic_contribution_vector <- t(genotype_matrix_causal) %*% beta_matrix
  
  # Compute Var(XB)
  var = var(genetic_contribution_vector)
  
  # Draw error values 
  error_vector <- rnorm(sample_size, mean = 0, sd = sqrt(1-h2))
  
  # Compute phenotype
  phenotype_vector <- genetic_contribution_vector + error_vector
  
  # Print relevant statistics for phenotype simulation 
  print(paste("Heritability:", var/var(phenotype_vector), sep=" "))
  
  # Create error file 
  error_file <- cbind(plink_data$fam$fam, plink_data$fam$id, genetic_contribution_vector, error_vector)
  colnames(error_file) <- c("id_1","id_2", "XB", "Error")
  
  # Create phenotype file 
  phenotype_file <- cbind(plink_data$fam$fam, plink_data$fam$id, phenotype_vector)

  # Write phenotype file
  print("Writing phenotype file...")
  write.table(phenotype_file, file = paste(output_file_prefix, ".phen", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=FALSE) 
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".phen", sep = ""), sep=" "))
  
  # Write par file
  print("Writing par file...")
  write.table(par_file, file = paste(output_file_prefix, ".par", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=TRUE)   
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".par", sep = ""), sep=" "))

  # Write error file 
  print("Writing error file...")
  write.table(error_file, file = paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep=" "))
}

target_linear_pheno_sim <- function(plink_data, sd_data, source_causal_loci, source_par_file, genetic_correlation, h2, output_file_prefix) {

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
  
  print(paste("Number of causal SNPs:", dim(genotype_matrix_causal)[1], sep=" "))
  print(paste("Sample Size:", dim(genotype_matrix_causal)[2], sep=" "))
  
  # Get reference allele from bim file of plink_data 
  bim <- plink_data$bim
  causal_bim <- bim[bim$id %in% causal_genotype_quote, ]
  ref_alleles <- causal_bim$ref
  
  # Read in source betas file and extract beta values
  par_file <- read.table(source_par_file, sep = "\t", header=TRUE)
  source_betas <- as.numeric(par_file$Effect)
  
  print("Starting beta simulation for causal SNPs...")
  
  # Draw beta values 
  betas <- rnorm(num_of_causal_snps, mean = source_betas, sd=sqrt((1-genetic_correlation)*h2/num_of_causal_snps))
  
  # Compute allele frequencies 
  allele_frequencies <- rowSums(genotype_matrix_causal) / (2 * ncol(genotype_matrix_causal))
  
  # Create par_file with relevant information 
  par_file <- cbind(rownames(genotype_matrix_causal), ref_alleles, allele_frequencies, betas)
  colnames(par_file) <- c("QTL", "RefAllele", "Frequency", "Effect")
  
  print("Starting phenotype simulation for each individual...")
  
  # Compute XB for each individual
  beta_matrix <- matrix(betas, nrow = length(betas), ncol = 1)
  
  genetic_contribution_vector <- t(genotype_matrix_causal) %*% beta_matrix
  
  # Compute Var(XB)
  var = var(genetic_contribution_vector)
  
  # Draw error values 
  error_vector <- rnorm(sample_size, mean = 0, sd = sqrt(1-h2))
  
  # Compute phenotype
  phenotype_vector <- genetic_contribution_vector + error_vector
  
  # Print relevant statistics for phenotype simulation 
  print(paste("Heritability:", var/var(phenotype_vector), sep=" "))
  
  # Create error file 
  error_file <- cbind(plink_data$fam$fam, plink_data$fam$id, genetic_contribution_vector, error_vector)
  colnames(error_file) <- c("id_1","id_2", "XB", "Error")
  
  # Create phenotype file 
  phenotype_file <- cbind(plink_data$fam$fam, plink_data$fam$id, phenotype_vector)

  # Write phenotype file
  print("Writing phenotype file...")
  write.table(phenotype_file, file = paste(output_file_prefix, ".phen", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=FALSE) 
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".phen", sep = ""), sep=" "))
  
  # Write par file
  print("Writing par file...")
  write.table(par_file, file = paste(output_file_prefix, ".par", sep = ""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names=TRUE)   
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, ".par", sep = ""), sep=" "))

  # Write error file 
  print("Writing error file...")
  write.table(error_file, file = paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
  print(paste("Sucessfully wrote file:", paste(output_file_prefix, "_pheno_error.txt", sep = ""), sep=" "))
}
