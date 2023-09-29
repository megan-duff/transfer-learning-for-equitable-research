#!/bin/bash
#SBATCH --job-name=TSI_linear_pheno_sim_redo
#SBATCH --output=TSI_linear_pheno_sim_%a_redo.out
#SBATCH --error=TSI_linear_pheno_sim_%a_redo.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=23,61,153,235,244,246,255,75

sim=${SLURM_ARRAY_TASK_ID}

ancestry="TSI"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation_redo.R "$ancestry" "$sim"
