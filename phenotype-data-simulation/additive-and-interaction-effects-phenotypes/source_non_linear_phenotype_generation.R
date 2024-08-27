#!/usr/bin/env Rscript

library(genio) 

source("/home/duffme/Scripts/pheno_sim_scripts/new_non_linear_pheno_sim_functions.R") 

subset=1

file="subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k" 

print("Reading in plink data...") 

plink_data<-load_in_data(file) 

bim_data<-get_bim(plink_data) 

source("/home/duffme/Scripts/pheno_sim_scripts/new_non_linear_pheno_sim_functions.R") 

print("Starting to standardize data...") 

sd_data<-sd_data(plink_data, "WG_af.afreq") 

print("Start phenotype simulation...") 

source_non_linear_pheno_sim_XB_ZY(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="1000_causal_snps.txt",  
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes") 

print("Finished: e = 0.05 / num_causal_snps = 1k / num_interaction_terms = 100") 

source_non_linear_pheno_sim_XB_ZY(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="1000_causal_snps.txt",  
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes")

print("Finished: e = 0.1 / num_causal_snps = 1k / num_interaction_terms = 100") 

source_non_linear_pheno_sim_XB_ZY(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="10000_causal_snps.txt",  
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes") 

print("Finished: e = 0.05 / num_causal_snps = 10k / num_interaction_terms = 100") 

source_non_linear_pheno_sim_XB_ZY(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="10000_causal_snps.txt",  
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes") 

print("Finished: e = 0.1 / num_causal_snps = 10k / num_interaction_terms = 100") 
print("---------------------------------------------------------------------------") 

print(paste("Finished subset:", subset)) 

rm(list = ls()) 

for (subset in 2:6){ 
    print(paste("Start subset:", subset)) 
    file=paste("subset_", subset, "_whole_genome_GBR_1000GP_hm3_Phase3_20k", sep="") 
    print("Reading in plink data...") 
    source("/home/duffme/Scripts/pheno_sim_scripts/new_non_linear_pheno_sim_functions.R") 
    plink_data<-load_in_data(file) 
    bim_data<-get_bim(plink_data) 
    print("Starting to standardize data...") 
    sd_data<-sd_data(plink_data, "WG_af.afreq") 
    print("Start phenotype simulation...") 

    source_non_linear_pheno_sim_XB_ZY_subset(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="1000_causal_snps.txt",  
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          source_additive_betas_file="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes_additive_betas.txt", 
                          source_interaction_file="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes_interactions.txt", 
                          output_file_prefix= paste("subset_", subset, "_GBR_1k_0.5_0.05_non_linear_phenotypes", sep="")) 

    print("Finished: e = 0.05 / num_causal_snps = 1k / num_interaction_terms = 100") 
    
    source_non_linear_pheno_sim_XB_ZY_subset(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="1000_causal_snps.txt",  
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          source_additive_betas_file="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes_additive_betas.txt", 
                          source_interaction_file="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes_interactions.txt", 
                          output_file_prefix= paste("subset_", subset, "_GBR_1k_0.5_0.1_non_linear_phenotypes", sep="")) 
    
    print("Finished: e = 0.1 / num_causal_snps = 1k / num_interaction_terms = 100") 
    
    source_non_linear_pheno_sim_XB_ZY_subset(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="10000_causal_snps.txt",  
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          source_additive_betas_file="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes_additive_betas.txt", 
                          source_interaction_file="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes_interactions.txt", 
                          output_file_prefix= paste("subset_", subset, "_GBR_10k_0.5_0.05_non_linear_phenotypes", sep="")) 
    
    print("Finished: e = 0.05 / num_causal_snps = 10k / num_interaction_terms = 100")
    
    source_non_linear_pheno_sim_XB_ZY_subset(sd_data=sd_data,  
                          plink_data=plink_data, 
                          causal_loci="10000_causal_snps.txt",  
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          source_additive_betas_file="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes_additive_betas.txt", 
                          source_interaction_file="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes_interactions.txt", 
                          output_file_prefix= paste("subset_", subset, "_GBR_10k_0.5_0.1_non_linear_phenotypes", sep="")) 
    
    print("Finished: e = 0.1 / num_causal_snps = 10k / num_interaction_terms = 100")

    print(paste("Finished subset:", subset)) 

    rm(list = ls()) 
} 

print("---------------------------------------------------------------------------") 

source("/home/duffme/Scripts/pheno_sim_scripts/new_non_linear_pheno_sim_functions.R")

print("Scale Effects")

source_non_linear_pheno_sim_scale_effects( 
  subset_1_data="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes",  
  subset_2_data="subset_2_GBR_1k_0.5_0.05_non_linear_phenotypes", 
  subset_3_data="subset_3_GBR_1k_0.5_0.05_non_linear_phenotypes", 
  subset_4_data="subset_4_GBR_1k_0.5_0.05_non_linear_phenotypes", 
  subset_5_data="subset_5_GBR_1k_0.5_0.05_non_linear_phenotypes",
  subset_6_data="subset_6_GBR_1k_0.5_0.05_non_linear_phenotypes", 
  causal_loci="1000_causal_snps.txt",  
  e=0.05,  
  num_interaction_terms=100, 
  h2=0.5, 
  sample_size=120000, 
  output_file_prefix="GBR_1k_0.5_0.05_non_linear_phenotypes", 
  source_par_file="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes_additive_betas.txt",  
  source_interaction_file="subset_1_GBR_1k_0.5_0.05_non_linear_phenotypes_interactions.txt") 

print("Finished: e = 0.05 / num_causal_snps = 1k / num_interaction_terms = 100")

source_non_linear_pheno_sim_scale_effects( 
  subset_1_data="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes",  
  subset_2_data="subset_2_GBR_1k_0.5_0.1_non_linear_phenotypes", 
  subset_3_data="subset_3_GBR_1k_0.5_0.1_non_linear_phenotypes", 
  subset_4_data="subset_4_GBR_1k_0.5_0.1_non_linear_phenotypes", 
  subset_5_data="subset_5_GBR_1k_0.5_0.1_non_linear_phenotypes",
  subset_6_data="subset_6_GBR_1k_0.5_0.1_non_linear_phenotypes", 
  causal_loci="1000_causal_snps.txt",  
  e=0.1,  
  num_interaction_terms=100, 
  h2=0.5, 
  sample_size=120000, 
  output_file_prefix="GBR_1k_0.5_0.1_non_linear_phenotypes", 
  source_par_file="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes_additive_betas.txt",  
  source_interaction_file="subset_1_GBR_1k_0.5_0.1_non_linear_phenotypes_interactions.txt") 

print("Finished: e = 0.1 / num_causal_snps = 1k / num_interaction_terms = 100")

source_non_linear_pheno_sim_scale_effects( 
  subset_1_data="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes",  
  subset_2_data="subset_2_GBR_10k_0.5_0.05_non_linear_phenotypes", 
  subset_3_data="subset_3_GBR_10k_0.5_0.05_non_linear_phenotypes", 
  subset_4_data="subset_4_GBR_10k_0.5_0.05_non_linear_phenotypes", 
  subset_5_data="subset_5_GBR_10k_0.5_0.05_non_linear_phenotypes",
  subset_6_data="subset_6_GBR_10k_0.5_0.05_non_linear_phenotypes", 
  causal_loci="10000_causal_snps.txt",  
  e=0.05,  
  num_interaction_terms=100, 
  h2=0.5, 
  sample_size=120000, 
  output_file_prefix="GBR_10k_0.5_0.05_non_linear_phenotypes", 
  source_par_file="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes_additive_betas.txt",  
  source_interaction_file="subset_1_GBR_10k_0.5_0.05_non_linear_phenotypes_interactions.txt") 

print("Finished: e = 0.05 / num_causal_snps = 10k / num_interaction_terms = 100")

source_non_linear_pheno_sim_scale_effects( 
  subset_1_data="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes",  
  subset_2_data="subset_2_GBR_10k_0.5_0.1_non_linear_phenotypes", 
  subset_3_data="subset_3_GBR_10k_0.5_0.1_non_linear_phenotypes", 
  subset_4_data="subset_4_GBR_10k_0.5_0.1_non_linear_phenotypes", 
  subset_5_data="subset_5_GBR_10k_0.5_0.1_non_linear_phenotypes",
  subset_6_data="subset_6_GBR_10k_0.5_0.1_non_linear_phenotypes", 
  causal_loci="10000_causal_snps.txt",  
  e=0.1,  
  num_interaction_terms=100, 
  h2=0.5, 
  sample_size=120000, 
  output_file_prefix="GBR_10k_0.5_0.1_non_linear_phenotypes", 
  source_par_file="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes_additive_betas.txt",  
  source_interaction_file="subset_1_GBR_10k_0.5_0.1_non_linear_phenotypes_interactions.txt") 

print("Finished: e = 0.1 / num_causal_snps = 10k / num_interaction_terms = 100")

print("---------------------------------------------------------------------------") 

print("Compute Phenotypes") 

for (subset in 1:6){ 
  print(paste("Start subset:", subset)) 
  file=paste("subset_", subset, "_whole_genome_GBR_1000GP_hm3_Phase3_20k", sep="") 
  print("Reading in plink data...") 
  source("/home/duffme/Scripts/pheno_sim_scripts/new_linear_pheno_sim_functions.R") 
  plink_data<-load_in_data(file) 
  bim_data<-get_bim(plink_data) 
  print("Starting to standardize data...") 
  sd_data<-sd_data(plink_data, "WG_af.afreq") 
  print("Start phenotype simulation...") 
  
  source_non_linear_pheno_sim_scaled_pheno_subset(
    sd_data=sd_data,  
    plink_data=plink_data, 
    causal_loci="1000_causal_snps.txt", 
    subset_num=subset, 
    scaled_additive_betas_file="GBR_1k_0.5_0.05_non_linear_phenotypes_scaled_additive_betas.txt", 
    scaled_interaction_file="GBR_1k_0.5_0.05_non_linear_phenotypes_scaled_interaction_betas.txt", 
    scaled_error_file="GBR_1k_0.5_0.05_non_linear_phenotypes_scaled_error_vector.txt", 
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix=paste("subset_",subset,"_GBR_1k_0.5_0.05_non_linear_phenotypes",sep="")) 
  
  print("Finished: e = 0.05 / num_causal_snps = 1k / num_interaction_terms = 100")

  source_non_linear_pheno_sim_scaled_pheno_subset(
    sd_data=sd_data,  
    plink_data=plink_data, 
    causal_loci="1000_causal_snps.txt", 
    subset_num=subset, 
    scaled_additive_betas_file="GBR_1k_0.5_0.1_non_linear_phenotypes_scaled_additive_betas.txt", 
    scaled_interaction_file="GBR_1k_0.5_0.1_non_linear_phenotypes_scaled_interaction_betas.txt", 
    scaled_error_file="GBR_1k_0.5_0.1_non_linear_phenotypes_scaled_error_vector.txt", 
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix=paste("subset_",subset,"_GBR_1k_0.5_0.1_non_linear_phenotypes",sep="")) 
  
  print("Finished: e = 0.1 / num_causal_snps = 1k / num_interaction_terms = 100")
  
  source_non_linear_pheno_sim_scaled_pheno_subset(
    sd_data=sd_data,  
    plink_data=plink_data, 
    causal_loci="10000_causal_snps.txt", 
    subset_num=subset, 
    scaled_additive_betas_file="GBR_10k_0.5_0.05_non_linear_phenotypes_scaled_additive_betas.txt", 
    scaled_interaction_file="GBR_10k_0.5_0.05_non_linear_phenotypes_scaled_interaction_betas.txt", 
    scaled_error_file="GBR_10k_0.5_0.05_non_linear_phenotypes_scaled_error_vector.txt", 
                          e=0.05,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix=paste("subset_",subset,"_GBR_10k_0.5_0.05_non_linear_phenotypes",sep="")) 
  
  print("Finished: e = 0.05 / num_causal_snps = 10k / num_interaction_terms = 100")

  source_non_linear_pheno_sim_scaled_pheno_subset(
    sd_data=sd_data,  
    plink_data=plink_data, 
    causal_loci="10000_causal_snps.txt", 
    subset_num=subset, 
    scaled_additive_betas_file="GBR_10k_0.5_0.1_non_linear_phenotypes_scaled_additive_betas.txt", 
    scaled_interaction_file="GBR_10k_0.5_0.1_non_linear_phenotypes_scaled_interaction_betas.txt", 
    scaled_error_file="GBR_10k_0.5_0.1_non_linear_phenotypes_scaled_error_vector.txt", 
                          e=0.1,  
                          num_interaction_terms=100, 
                          h2=0.5,  
                          output_file_prefix=paste("subset_",subset,"_GBR_10k_0.5_0.1_non_linear_phenotypes",sep="")) 
  
  print("Finished: e = 0.1 / num_causal_snps = 10k / num_interaction_terms = 100")
  
}

print("---------------------------------------------------------------------------")

print(paste("Finished simulation!!!")) 
