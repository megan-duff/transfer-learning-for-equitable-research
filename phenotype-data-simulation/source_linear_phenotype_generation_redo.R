#!/usr/bin/env Rscript

library(genio)

source("/home/duffme/Scripts/pheno_sim_scripts/new_linear_pheno_sim_functions.R")

args <- commandArgs(trailingOnly = TRUE)

sim <- as.integer(args[1])

subset=1

file_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_", subset, "_GBR_100k_0.8_phenotypes.phen", sep="")

# Check if the file exists
if (file.exists(file_path)) {
  print(paste("Subset", subset, "file exists. No need to redo...", sep=""))
} else {
  print(paste("Start subset:", subset))
  
  file=paste("subset_", subset, "_whole_genome_GBR_1000GP_hm3_Phase3_20k", sep="")
  
  print("Reading in plink data...")
  plink_data<-load_in_data(file)
  
  bim_data<-get_bim(plink_data)
  
  source("/home/duffme/Scripts/pheno_sim_scripts/new_linear_pheno_sim_functions.R")
  
  print("Starting to standardize data...")
  sd_data<-sd_data(plink_data, "WG_af.afreq")
  
  print("Start phenotype simulation...")

  source_linear_pheno_sim(sd_data=sd_data, 
                          plink_data=plink_data,
                          causal_loci="100_causal_snps.txt", 
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 100")

  source_linear_pheno_sim(sd_data=sd_data, 
                          plink_data=plink_data,
                          causal_loci="1000_causal_snps.txt", 
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_1k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 1k")

  source_linear_pheno_sim(sd_data=sd_data, 
                          plink_data=plink_data,
                          causal_loci="10000_causal_snps.txt", 
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_10k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 10k")

  source_linear_pheno_sim(sd_data=sd_data,
                          plink_data=plink_data,
                          causal_loci="1e+05_causal_snps.txt", 
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 100k")

  source_linear_pheno_sim(sd_data=sd_data,
                          plink_data=plink_data,
                          causal_loci="100_causal_snps.txt", 
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 100")

  source_linear_pheno_sim(sd_data=sd_data,
                          plink_data=plink_data,
                          causal_loci="1000_causal_snps.txt", 
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_1k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 1k")

  source_linear_pheno_sim(sd_data=sd_data,
                          plink_data=plink_data,
                          causal_loci="10000_causal_snps.txt", 
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_10k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 10k")

  source_linear_pheno_sim(sd_data=sd_data,
                          plink_data=plink_data,
                          causal_loci="1e+05_causal_snps.txt", 
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 100k")
  
  print(paste("Finished subset:", subset))

  rm(list = ls())
  }

for (subset in 2:6){

  file_path = paste("/scratch/duffme_gensim/Simulations/GBR/sim_", sim, "/subset_", subset, "_GBR_100k_0.8_phenotypes.phen", sep="")

  # Check if the file exists
  if (file.exists(file_path)) {
  print(paste("Subset", subset, "file exists. No need to redo...", sep=""))
} else {
  print(paste("Start subset:", subset))
  
  file=paste("subset_", subset, "_whole_genome_GBR_1000GP_hm3_Phase3_20k", sep="")
  
  print("Reading in plink data...")
  
  source("/home/duffme/Scripts/pheno_sim_scripts/new_linear_pheno_sim_functions.R")
  
  plink_data<-load_in_data(file)
  
  bim_data<-get_bim(plink_data)
  
  print("Starting to standardize data...")
  sd_data<-sd_data(plink_data, "WG_af.afreq")
  
  print("Start phenotype simulation...")

  source_linear_pheno_sim_subset(sd_data=sd_data, 
                          plink_data=plink_data,
                          source_causal_loci="100_causal_snps.txt", 
                          source_par_file="subset_1_GBR_100_0.4_phenotypes.par",
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 100")

  source_linear_pheno_sim_subset(sd_data=sd_data, 
                          plink_data=plink_data,
                          source_causal_loci="1000_causal_snps.txt", 
                          source_par_file="subset_1_GBR_1k_0.4_phenotypes.par",
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_1k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 1k")

  source_linear_pheno_sim_subset(sd_data=sd_data, 
                          plink_data=plink_data,
                          source_causal_loci="10000_causal_snps.txt", 
                          source_par_file="subset_1_GBR_10k_0.4_phenotypes.par",
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_10k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 10k")

  source_linear_pheno_sim_subset(sd_data=sd_data,
                          plink_data=plink_data,
                          source_causal_loci="1e+05_causal_snps.txt", 
                          source_par_file="subset_1_GBR_100k_0.4_phenotypes.par",
                          h2=0.4, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100k_0.4_phenotypes", sep=""))
  
  print("Finished: h2 = 0.4 and num_causal_snps = 100k")

  source_linear_pheno_sim_subset(sd_data=sd_data,
                          plink_data=plink_data,
                          source_causal_loci="100_causal_snps.txt",
                          source_par_file="subset_1_GBR_100_0.8_phenotypes.par",
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 100")

  source_linear_pheno_sim_subset(sd_data=sd_data,
                          plink_data=plink_data,
                          source_causal_loci="1000_causal_snps.txt",
                          source_par_file="subset_1_GBR_1k_0.8_phenotypes.par",
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_1k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 1k")

  source_linear_pheno_sim_subset(sd_data=sd_data,
                          plink_data=plink_data,
                          source_causal_loci="10000_causal_snps.txt",
                          source_par_file="subset_1_GBR_10k_0.8_phenotypes.par",
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_10k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 10k")

  source_linear_pheno_sim_subset(sd_data=sd_data,
                          plink_data=plink_data,
                          source_causal_loci="1e+05_causal_snps.txt",
                          source_par_file="subset_1_GBR_100k_0.8_phenotypes.par",
                          h2=0.8, 
                          output_file_prefix=paste("subset_", subset, "_GBR_100k_0.8_phenotypes", sep=""))
  
  print("Finished: h2 = 0.8 and num_causal_snps = 100k")
  
  print(paste("Finished subset:", subset))
  
  rm(list = ls())
  }
}

print("---------------------------------------------------------------------------")
print(paste("Finished simulation redo!!!"))

