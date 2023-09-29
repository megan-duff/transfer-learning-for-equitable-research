#!/bin/bash
#SBATCH --job-name=CEU_linear_pheno_sim
#SBATCH --output=CEU_linear_pheno_sim_%a.out
#SBATCH --error=CEU_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

ancestry="CEU"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation.R "$ancestry" "$sim"
