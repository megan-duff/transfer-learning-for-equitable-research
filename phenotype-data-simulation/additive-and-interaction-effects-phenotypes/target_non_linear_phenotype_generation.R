#!/usr/bin/env Rscript

library(genio)

source("/home/duffme/Scripts/pheno_sim_scripts/new_non_linear_pheno_sim_functions.R")

# Grab ancestry and sim variables passed from bash script that are needed to run following code
args <- commandArgs(trailingOnly = TRUE)
ancestry <- args[1]
sim <- args[2]

# Define whole genome full sample data set name
whole_data_file = paste("whole_genome_", ancestry, "_1000GP_hm3_Phase3_20k", sep="")

# Read in bed/bim/fam plink files for target data set 
plink_data<-load_in_data(whole_data_file)

# Standardize genotype data 
sd_data <- sd_data(plink_data, "WG_af.afreq")

# Define paths are various files called during phenotype simulation algorithm 
source_causal_loci_1k_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/1000_causal_snps.txt", sep="")
source_causal_loci_10k_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/10000_causal_snps.txt", sep="")

source_par_file_1k_0.05_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes.par", sep="")
source_par_file_1k_0.1_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes.par", sep="")

source_par_file_10k_0.05_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes.par", sep="")
source_par_file_10k_0.1_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes.par", sep="")

source_interaction_file_1k_0.05_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes_interactions.txt", sep="")
source_interaction_file_1k_0.1_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes_interactions.txt", sep="")

source_interaction_file_10k_0.05_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes_interactions.txt", sep="")
source_interaction_file_10k_0.1_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes_interactions.txt", sep="")

# Simulate phenotype 1: h2=0.5, num_causal_snps=1k, num_interaction_terms=100, e=0.05
target_non_linear_pheno_sim(sd_data=sd_data,
                        plink_data=plink_data,
                        source_causal_loci=source_causal_loci_1k_path, 
                        source_par_file=source_par_file_1k_0.05_path, 
                        source_interaction_file=source_interaction_file_1k_0.05_path,
                        genetic_correlation=0.8, 
                        h2=0.5, 
                        e=0.05,
                        num_interaction_terms=100,
                        output_file_prefix=paste(ancestry,"_1k_0.5_0.05_non_linear_phenotypes",sep=""))

# Simulate phenotype 2: h2=0.5, num_causal_snps=1k, num_interaction_terms=100, e=0.1
target_non_linear_pheno_sim(sd_data=sd_data,
                        plink_data=plink_data,
                        source_causal_loci=source_causal_loci_1k_path, 
                        source_par_file=source_par_file_1k_0.1_path, 
                        source_interaction_file=source_interaction_file_1k_0.1_path,
                        genetic_correlation=0.8, 
                        h2=0.5, 
                        e=0.1,
                        num_interaction_terms=100,
                        output_file_prefix=paste(ancestry,"_1k_0.5_0.1_non_linear_phenotypes",sep=""))

# Simulate phenotype 3: h2=0.5, num_causal_snps=10k, num_interaction_terms=100, e=0.05
target_non_linear_pheno_sim(sd_data=sd_data,
                        plink_data=plink_data,
                        source_causal_loci=source_causal_loci_10k_path, 
                        source_par_file=source_par_file_10k_0.05_path, 
                        source_interaction_file=source_interaction_file_10k_0.05_path,
                        genetic_correlation=0.8, 
                        h2=0.5, 
                        e=0.05,
                        num_interaction_terms=100,
                        output_file_prefix=paste(ancestry,"_10k_0.5_0.05_non_linear_phenotypes",sep=""))

# Simulate phenotype 4: h2=0.5, num_causal_snps=10k, num_interaction_terms=100, e=0.1
target_non_linear_pheno_sim(sd_data=sd_data,
                        plink_data=plink_data,
                        source_causal_loci=source_causal_loci_10k_path, 
                        source_par_file=source_par_file_10k_0.1_path, 
                        source_interaction_file=source_interaction_file_10k_0.1_path,
                        genetic_correlation=0.8, 
                        h2=0.5, 
                        e=0.1,
                        num_interaction_terms=100,
                        output_file_prefix=paste(ancestry,"_10k_0.5_0.1_non_linear_phenotypes",sep=""))
