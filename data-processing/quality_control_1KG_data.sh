#!/bin/bash

#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out WG_1000G_phase3

#Keep only autosomal SNPs, extract hapmap3 variants, remove first and second degree relatives, and maf filter of 0.05 for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --maf 0.05 --out WG_1000G_phase3_common

#Keep only autosomal SNPs, extract hapmap3 variants, remove first and second degree relatives, maf filter of 0.05, and LD prune for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --indep-pairwise 50 5 0.5 --maf 0.05 --out WG_1000G_phase3_HM3_common_ldpruned

#rename .psam file to have the same prefix as the .pgen files (this makes it easier to input into plink)
for chr in {1..22}; do cp all_phase3.psam chr_${chr}.psam; done

#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for chromosome specific files
for chr in {1..22}; do 
../.././plink2 --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do mv chr_${chr}.psam chr${chr}_1000G_phase3.psam ; done

#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for chromosome specific files & save the file as .haps/.legend for HAPGEN2
for chr in {1..22}; do 
../.././plink2 --export hapslegend --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out chr${chr}_1000G_phase3_HM3; done
