#!/bin/bash
#SBATCH --job-name=GBR_snp_list
#SBATCH --output=GBR_snp_list_%a.out
#SBATCH --error=GBR_snp_list_%a.err
#SBATCH -n 5
#SBATCH -p math-alderaan
#SBATCH --array 1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/causal_snp_generation.R ${sim}
