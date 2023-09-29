#!/bin/bash
#SBATCH --job-name=GBR_train_val_test_file_generation_phenotype_files_only
#SBATCH --output=GBR_linear_pheno_sim_%a_phenotype_files_only.out
#SBATCH --error=GBR_linear_pheno_sim_%a_phenotype_files_only.err
#SBATCH --partition=math-alderaan
#SBATCH -n 5
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/source_train_val_test_phenotype_files.R
