#!/bin/bash
#SBATCH --job-name=CHB_af
#SBATCH --output=CHB_af
#SBATCH --error=CHB_af_err
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/CHB/sim_${sim}

for chr in {1..22}; do ../../.././plink2 --bfile CHB_1000GP_hm3_Phase3_chr${chr}_20k.controls --freq --out chr_${chr}_af; done

cat chr_*_af.afreq > sim_${sim}_WG_af.txt
