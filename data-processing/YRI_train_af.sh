#!/bin/bash
#SBATCH --job-name=YRI_train_af
#SBATCH --output=YRI_train_af_%a.out
#SBATCH --error=YRI_train_af_%a.err
#SBATCH -n 5
#SBATCH -p math-alderaan
#SBATCH --array=1-100

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/YRI/sim_${sim}

/home/duffme/Softwares/plink2 --bfile whole_genome_YRI_1000GP_hm3_Phase3_20k_train --freq --out train_WG_af

echo "Done with sim ${sim}"
