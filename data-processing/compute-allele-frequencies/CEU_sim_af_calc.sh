#!/bin/bash
#SBATCH --job-name=CEU_af
#SBATCH --output=CEU_af
#SBATCH --error=CEU_af_err
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/CEU/sim_${sim}

for chr in {1..22}; do ../../.././plink2 --bfile CEU_1000GP_hm3_Phase3_chr${chr}_20k.controls --freq --out chr_${chr}_af; done

cat chr_*_af.afreq > sim_${sim}_WG_af.txt
