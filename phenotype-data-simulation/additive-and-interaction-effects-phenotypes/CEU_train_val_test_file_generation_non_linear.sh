#!/bin/bash
#SBATCH --job-name=CEU_train_val_test_file_generation_non_linear
#SBATCH --output=CEU_non_linear_pheno_sim_%a.out
#SBATCH --error=CEU_non_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -n 10
#SBATCH --array=1-110

ancestry="CEU"

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_train_val_test_phenotype_files_non_linear.R "$ancestry"
