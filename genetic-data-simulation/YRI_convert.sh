#!/bin/bash
#SBATCH --job-name=YRI_convert
#SBATCH --output=YRI_convert
#SBATCH --error=YRI_convert
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/YRI/sim_${sim}

for chr in {1..22};
  do
  ../../.././plink --data ASW_1000GP_hm3_Phase3_chr${chr}_20k.controls\
		--make-bed\
		--allow-no-sex\
		--oxford-single-chr ${chr}\
		--out YRI_1000GP_hm3_Phase3_chr${chr}_20k.controls;
  done
