#!/bin/bash
#SBATCH --job-name=scaling_train_val_test_data
#SBATCH --output=scaling_train_val_test_data.out
#SBATCH --error=scaling_train_val_test_data.err
#SBATCH -n 64
#SBATCH -p math-alderaan
#SBATCH --array=1

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}/

singularity exec /storage/singularity/tensorflow.sif ~/Scripts/pheno_sim_scripts/scaling_train_val_test_data.py
