#!/bin/bash
#SBATCH --job-name=GBR_train_af
#SBATCH --output=GBR_train_af_%a.out
#SBATCH --error=GBR_train_af_%a.err
#SBATCH -n 5
#SBATCH -p math-alderaan
#SBATCH --array=11-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

../../.././plink2 --bfile whole_genome_GBR_1000GP_hm3_Phase3_120k_train --freq --out train_WG_af

echo "Done with sim ${sim}"
