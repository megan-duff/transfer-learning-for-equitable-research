#!/bin/bash
#SBATCH --job-name=CEU_train_val_test_file_generation
#SBATCH --output=CEU_linear_pheno_sim_%a.out
#SBATCH --error=CEU_linear_pheno_sim_%a.err
#SBATCH --partition=math-alderaan
#SBATCH -n 10
#SBATCH --array=1-100

ancestry="CEU"

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}/

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_train_val_test_id_files.R 

singularity exec /storage/singularity/mixtures.sif ~/Scripts/pheno_sim_scripts/target_train_val_test_phenotype_files.R "$ancestry"

/home/duffme/Softwares/plink --bfile whole_genome_${ancestry}_1000GP_hm3_Phase3_20k --keep train_samples.txt --make-bed --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k_train

/home/duffme/Softwares/plink --bfile whole_genome_${ancestry}_1000GP_hm3_Phase3_20k --keep val_samples.txt --make-bed --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k_val

/home/duffme/Softwares/plink --bfile whole_genome_${ancestry}_1000GP_hm3_Phase3_20k --keep test_samples.txt --make-bed --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k_test
