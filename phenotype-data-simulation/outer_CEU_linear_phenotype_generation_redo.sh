#!/bin/bash
#SBATCH --job-name=CEU_linear_pheno_sim_redo
#SBATCH --output=CEU_linear_pheno_sim_%a_redo.out
#SBATCH --error=CEU_linear_pheno_sim_%a_redo.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=35,56,57,132,169,176,224,339,150,336,347

sim=${SLURM_ARRAY_TASK_ID}

ancestry="CEU"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation_redo.R "$ancestry" "$sim"
