#!/bin/bash
#SBATCH --job-name=GBR_non_linear_pheno_sim
#SBATCH --output=GBR_non_linear_pheno_sim_%a.out
#SBATCH --error=GBR_non_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --array=10-110%10

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/source_non_linear_phenotype_generation.R
