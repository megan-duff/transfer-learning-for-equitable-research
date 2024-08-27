#!/bin/bash
#SBATCH --job-name=CHB_linear_pheno_sim_redo
#SBATCH --output=CHB_linear_pheno_sim_%a_redo.out
#SBATCH --error=CHB_linear_pheno_sim_%a_redo.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=1,15,30,36,39,75,100,105,136,137,150,157,160,172,194,325

sim=${SLURM_ARRAY_TASK_ID}

ancestry="CHB"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_linear_phenotype_generation_redo.R "$ancestry" "$sim"
