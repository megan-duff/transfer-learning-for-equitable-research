#!/bin/bash
#SBATCH --job-name=YRI_gw_files
#SBATCH --output=YRI_gw_files
#SBATCH --error=YRI_gw_files_err
#SBATCH -n 1
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

ancestry="YRI"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}

for i in {2..22}; do
    echo "${ancestry}_1000GP_hm3_Phase3_chr${i}_20k.controls" >> merge-list.txt
done

plink --bfile ${ancestry}_1000GP_hm3_Phase3_chr1_20k.controls --merge-list merge-list.txt --make-bed --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k

for chr in {1..22}; do plink --bfile ${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls\
                --exclude whole_genome_${ancestry}_1000GP_hm3_Phase3_20k-merge.missnp\
                --make-bed \
                --out ${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k_missnp.controls; done

for chr in {1..22}; do echo ${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k_missnp.controls >> chr_nodup_file_list.txt; done

../../.././plink --merge-list chr_nodup_file_list.txt --allow-no-sex --make-bed --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k
