#!/usr/bin/env Rscript

library(genio)

source("/home/duffme/Scripts/pheno_sim_scripts/new_linear_pheno_sim_functions.R")

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
source_causal_loci_100_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/100_causal_snps.txt", sep="")
source_causal_loci_1k_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/1000_causal_snps.txt", sep="")
source_causal_loci_10k_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/10000_causal_snps.txt", sep="")
source_causal_loci_100k_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/1e+05_causal_snps.txt", sep="")

source_par_file_100_0.4_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_100_0.4_phenotypes.par", sep="")
source_par_file_100_0.8_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_100_0.8_phenotypes.par", sep="")

source_par_file_1k_0.4_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.4_phenotypes.par", sep="")
source_par_file_1k_0.8_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_1k_0.8_phenotypes.par", sep="")

source_par_file_10k_0.4_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.4_phenotypes.par", sep="")
source_par_file_10k_0.8_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_10k_0.8_phenotypes.par", sep="")

source_par_file_100k_0.4_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_100k_0.4_phenotypes.par", sep="")
source_par_file_100k_0.8_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_1_GBR_100k_0.8_phenotypes.par", sep="")

# Simulate phenotype 1: h2=0.4, num_causal_snps=100
file_path <- paste(ancestry, "_100_0.4_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data = sd_data,
    plink_data = plink_data,
    source_causal_loci = source_causal_loci_100_path,
    source_par_file = source_par_file_100_0.4_path,
    genetic_correlation = 0.8,
    h2 = 0.4,
    output_file_prefix = paste(ancestry, "_100_0.4_phenotypes", sep = ""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}

# Simulate phenotype 2: h2=0.4, num_causal_snps=1k
file_path <- paste(ancestry, "_1k_0.4_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data, 
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_1k_path, 
    source_par_file=source_par_file_1k_0.4_path, 
    genetic_correlation=0.8, 
    h2=0.4, 
    output_file_prefix=paste(ancestry,"_1k_0.4_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}


# Simulate phenotype 3: h2=0.4, num_causal_snps=10k
file_path <- paste(ancestry, "_10k_0.4_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data,
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_10k_path, 
    source_par_file=source_par_file_10k_0.4_path, 
    genetic_correlation=0.8, 
    h2=0.4, 
    output_file_prefix=paste(ancestry,"_10k_0.4_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}


# Simulate phenotype 4: h2=0.4, num_causal_snps=100k
file_path <- paste(ancestry, "_100k_0.4_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data, 
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_100k_path, 
    source_par_file=source_par_file_100k_0.4_path, 
    genetic_correlation=0.8, 
    h2=0.4, 
    output_file_prefix=paste(ancestry,"_100k_0.4_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}


# Simulate phenotype 5: h2=0.8, num_causal_snps=100
file_path <- paste(ancestry, "_100_0.8_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data,
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_100_path, 
    source_par_file=source_par_file_100_0.8_path, 
    genetic_correlation=0.8, 
    h2=0.8, 
    output_file_prefix=paste(ancestry,"_100_0.8_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}


# Simulate phenotype 6: h2=0.8, num_causal_snps=1k
file_path <- paste(ancestry, "_1k_0.8_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data,
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_1k_path, 
    source_par_file=source_par_file_1k_0.8_path, 
    genetic_correlation=0.8, 
    h2=0.8, 
    output_file_prefix=paste(ancestry,"_1k_0.8_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}

# Simulate phenotype 7: h2=0.8, num_causal_snps=10k
file_path <- paste(ancestry, "_10k_0.8_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data,
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_10k_path, 
    source_par_file=source_par_file_10k_0.8_path, 
    genetic_correlation=0.8, 
    h2=0.8, 
    output_file_prefix=paste(ancestry,"_10k_0.8_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}


# Simulate phenotype 8: h2=0.8, num_causal_snps=100k
file_path <- paste(ancestry, "_100k_0.8_phenotypes.phen", sep = "")
if (!file.exists(file_path)) {
  print("File not found & starting redo!")
  target_linear_pheno_sim(
    sd_data=sd_data,
    plink_data=plink_data,
    source_causal_loci=source_causal_loci_100k_path,
    source_par_file=source_par_file_100k_0.8_path, 
    genetic_correlation=0.8, 
    h2=0.8, 
    output_file_prefix=paste(ancestry,"_100k_0.8_phenotypes",sep=""))
  print(paste("File finished:", file_path,sep=" "))
} else {print(paste("File already existed:", file_path,sep=" "))}
