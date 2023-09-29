#!/bin/bash
#SBATCH --job-name=GBR_linear_pheno_sim_redo
#SBATCH --output=GBR_linear_pheno_sim_%a.out
#SBATCH --error=GBR_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan-gpu
#SBATCH -n 10
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/source_linear_phenotype_generation_redo.R "$sim"
