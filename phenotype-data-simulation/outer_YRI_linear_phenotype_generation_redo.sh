#!/bin/bash
#SBATCH --job-name=YRI_linear_pheno_sim_redo
#SBATCH --output=YRI_linear_pheno_sim_%a_redo.out
#SBATCH --error=YRI_linear_pheno_sim_%a_redo.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=3,8,66,83,257,261,370,383

sim=${SLURM_ARRAY_TASK_ID}

ancestry="YRI"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation_redo.R "$ancestry" "$sim"
