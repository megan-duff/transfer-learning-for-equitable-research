#!/bin/bash
#SBATCH --job-name=CHB_linear_pheno_sim
#SBATCH --output=CHB_linear_pheno_sim_%a.out
#SBATCH --error=CHB_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=75,98,325

sim=${SLURM_ARRAY_TASK_ID}

ancestry="CHB"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation.R "$ancestry" "$sim"
