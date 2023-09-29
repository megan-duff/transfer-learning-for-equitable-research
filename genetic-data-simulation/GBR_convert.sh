#!/bin/bash
#SBATCH --job-name=GBR_convert
#SBATCH --output=GBR_convert
#SBATCH --error=GBR_convert
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

for chr in {1..22};
  do
  ../../.././plink --data GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls\
		--make-bed\
		--allow-no-sex\
		--oxford-single-chr ${chr}\
		--out GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls;
  done
