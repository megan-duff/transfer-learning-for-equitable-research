#!/usr/bin/env Rscript

snp_list = c("100", "1k", "10k", "100k")
h2_list = c("0.4", "0.8")

train_samples<-read.table('train_samples.txt', header=TRUE)
val_samples<-read.table('val_samples.txt', header=TRUE)
test_samples<-read.table('test_samples.txt', header=TRUE)

for (snp in snp_list) {
  for (h2 in h2_list) {
    s_1=read.table(paste("subset_1_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")
    s_2=read.table(paste("subset_2_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")
    s_3=read.table(paste("subset_3_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")
    s_4=read.table(paste("subset_4_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")
    s_5=read.table(paste("subset_5_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")
    s_6=read.table(paste("subset_6_GBR_", snp,"_", h2, "_phenotypes.phen",sep=""), header = F, sep="\t")

    pheno = rbind(s_1, s_2, s_3, s_4, s_5, s_6)
    
    colnames(pheno) = c("ID", "ID2", "Phenotype")
    
    pheno_train <- pheno[pheno$ID %in% train_samples$V1,]
    pheno_val <- pheno[pheno$ID %in% val_samples$V1,]
    pheno_test <- pheno[pheno$ID %in% test_samples$V1,]
    
    pheno_train_file = paste("GBR_train_", snp, "_", h2, "_phenotypes.phen",sep="")
    pheno_val_file = paste("GBR_val_", snp, "_", h2, "_phenotypes.phen",sep="")
    pheno_test_file = paste("GBR_test_", snp, "_", h2, "_phenotypes.phen",sep="")
    
    write.table(pheno_train, file = pheno_train_file, row.names = FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    write.table(pheno_val, file = pheno_val_file, row.names = FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    write.table(pheno_test, file = pheno_test_file, row.names = FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    
  }
}
