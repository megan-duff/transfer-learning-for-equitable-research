#!/bin/bash
#SBATCH --job-name=clean_up_chr_files
#SBATCH --output=clean_up_chr_files_%a.out
#SBATCH --error=clean_up_chr_files_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -n 5
#SBATCH --array=1-22

chr=${SLURM_ARRAY_TASK_ID}

echo "Start chromosome ${chr}..."

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls.bim" -exec rm {} +

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.bim" -exec rm {} +

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.bed" -exec rm {} +

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls.bed" -exec rm {} +

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls.fam" -exec rm {} +

find /scratch/duffme_gensim/Simulations/GBR/ -type f -name "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.fam" -exec rm {} +

echo "Finished chromosome ${chr}!"
