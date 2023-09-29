#!/bin/bash
#SBATCH --job-name=TSI_clean_up_chr_files
#SBATCH --output=TSI_clean_up_chr_files_%a.out
#SBATCH --error=TSI_clean_up_chr_files_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -n 3
#SBATCH --array=1-22

ancestry="TSI"

chr=${SLURM_ARRAY_TASK_ID}

echo "Start chromosome ${chr}..."

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls.bim" -exec rm {} +

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls_no_missnp.bim" -exec rm {} +

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls_no_missnp.bed" -exec rm {} +

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls.bed" -exec rm {} +

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls.fam" -exec rm {} +

find /scratch/duffme_gensim/Simulations/${ancestry}/ -type f -name "${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls_no_missnp.fam" -exec rm {} +

echo "Finished chromosome ${chr}!"
