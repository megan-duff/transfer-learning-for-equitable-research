#!/bin/bash
#SBATCH --job-name=GBR_ref_af
#SBATCH --output=GBR_ref_af
#SBATCH --error=GBR_ref_af_err
#SBATCH -n 1
#SBATCH -p math-alderaan

sim=${SLURM_ARRAY_TASK_ID}

cd /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/

for chr in {1..22}; do ../../.././plink2 --bfile GBR_chr${chr}_1000G_phase3_HM3 --freq --out chr_${chr}_af; done

cat chr_*_af.afreq > ref_WG_af.txt
