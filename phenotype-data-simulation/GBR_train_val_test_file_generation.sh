#!/bin/bash

#SBATCH --job-name=GBR_train_val_test_file_generation
#SBATCH --output=GBR_linear_pheno_sim_%a.out
#SBATCH --error=GBR_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 5
#SBATCH --array=1-10

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/source_train_val_test_id_files.R

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/source_train_val_test_phenotype_files.R

../.././plink  --bfile whole_genome_GBR_1000GP_hm3_Phase3_120k --keep train_samples.txt --make-bed --out whole_genome_GBR_1000GP_hm3_Phase3_120k_train

../.././plink  --bfile whole_genome_GBR_1000GP_hm3_Phase3_120k --keep val_samples.txt --make-bed --out whole_genome_GBR_1000GP_hm3_Phase3_120k_val

../.././plink  --bfile whole_genome_GBR_1000GP_hm3_Phase3_120k --keep test_samples.txt --make-bed --out whole_genome_GBR_1000GP_hm3_Phase3_120k_test

if [[ -e whole_genome_GBR_1000GP_hm3_Phase3_120k_train.bed && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_train.bim && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_train.fam && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_val.bed && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_val.bim && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_val.fam && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_test.bed && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_test.bim && -e whole_genome_GBR_1000GP_hm3_Phase3_120k_test.fam ]]; then 
  echo "All files created"
  rm whole_genome_GBR_1000GP_hm3_Phase3_120k.bed
  rm whole_genome_GBR_1000GP_hm3_Phase3_120k.bim
  rm whole_genome_GBR_1000GP_hm3_Phase3_120k.fam
  echo "Successfully removed files - DONE"
else
  echo "Files missing -- need to redo!!"
fi
