#!/bin/bash

#SBATCH --job-name=pc_glm
#SBATCH --output=pc_glm_%a.out
#SBATCH --error=pc_glm_%a.err
#SBATCH -n 10
#SBATCH -p math-alderaan-gpu
#SBATCH --array=11-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

train_file="whole_genome_GBR_1000GP_hm3_Phase3_120k_train"

# Compute PCs for the training sample's genotypes
../../.././plink2 --bfile ${train_file} --pca approx --out ${train_file}

for h2 in 0.4 0.8; do
  for snp in "100" "1k" "10k" "100k"; do
    #Start phenotype GWAS
    echo "For h2: $h2 and snp: $snp"

    phenotype_file="GBR_train_${snp}_${h2}_phenotypes.phen"
    output_file="GBR_train_${snp}_${h2}_phenotypes"

    echo "FID IID PHENO1" | cat - ${phenotype_file} > temp && mv temp ${phenotype_file}

    # Run GWAS analysis for genotypes and corresponding phenotype
    ../../.././plink2 --bfile ${train_file} \
            --pheno ${phenotype_file} \
            --glm --out ${output_file} \
            --covar ${train_file}.eigenvec \
            --pheno-name PHENO1

    sort -k 9,9n ${output_file}.PHENO1.glm.linear > ${output_file}.sorted

    head -n 10001 ${output_file}.sorted | tail -n 10000 > ${output_file}_top_10k_snps.txt
    head -n 50001 ${output_file}.sorted | tail -n 50000 > ${output_file}_top_50k_snps.txt
    head -n 100001 ${output_file}.sorted | tail -n 100000 > ${output_file}_top_100k_snps.txt
    head -n 250001 ${output_file}.sorted | tail -n 250000 > ${output_file}_top_250k_snps.txt

  done
done
