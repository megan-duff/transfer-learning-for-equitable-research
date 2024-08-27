#!/bin/bash
#SBATCH --job-name=YRI_non_linear_pheno_sim
#SBATCH --output=YRI_non_linear_pheno_sim_%a.out
#SBATCH --error=YRI_non_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --array=1-110%10

sim=${SLURM_ARRAY_TASK_ID}

ancestry="YRI"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_non_linear_phenotype_generation.R "$ancestry" "$sim"
