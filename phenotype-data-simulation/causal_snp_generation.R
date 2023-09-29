#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

sim <- as.integer(args[1])

GBR_af <- read.table("WG_af.afreq")
TSI_af <- read.table(paste("/scratch/duffme_gensim/Simulations/TSI/sim_", sim, "/WG_af.afreq", sep=""))
YRI_af <- read.table(paste("/scratch/duffme_gensim/Simulations/YRI/sim_", sim, "/WG_af.afreq", sep=""))
CHB_af <- read.table(paste("/scratch/duffme_gensim/Simulations/CHB/sim_", sim, "/WG_af.afreq", sep=""))
CEU_af <- read.table(paste("/scratch/duffme_gensim/Simulations/CEU/sim_", sim, "/WG_af.afreq", sep=""))

master_af_df <- as.data.frame(cbind(GBR_af[,2], GBR_af[,5], TSI_af[,5], YRI_af[,5], CHB_af[,5], CEU_af[,5]))

colnames(master_af_df) <- c("rsID", "GBR_af", "TSI_af", "YRI_af", "CHB_af", "CEU_af")

master_af_df_remove_zero <- subset(master_af_df, !apply(master_af_df[, 2:6], 1, function(x) any(x == 0)))

snp_list = c(100, 1000, 10000, 100000)

for (num_of_snps in snp_list){
  sample = sample(master_af_df_remove_zero[,1], num_of_snps)
  write.table(sample, file=paste(num_of_snps,"_causal_snps.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

print("Finished causal SNP selection!")
