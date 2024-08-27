#!/bin/bash
#SBATCH --job-name=CHB_gw_af
#SBATCH --output=CHB_gw_af_%a.out
#SBATCH --error=CHB_gw_af_%a.err
#SBATCH -n 5
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

ancestry="CHB"

cd /scratch/duffme_gensim/Simulations/${ancestry}/sim_${sim}

for chr in {1..22}; do 
../../plink --bfile ${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls \
--extract /scratch/duffme_gensim/Simulations/GBR/sim_1/whole_genome_GBR_1000GP_hm3_Phase3_120k_train.bim \
--make-bed \
--allow-no-sex \
--out pre_merge_${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls;
done 

rm merge-list.txt

for chr in {2..22}; do echo "pre_merge_${ancestry}_1000GP_hm3_Phase3_chr${chr}_20k.controls" >> merge-list.txt; done

if [ ! -s "whole_genome_${ancestry}_1000GP_hm3_Phase3_20k.bed" ] || [ ! -s "whole_genome_${ancestry}_1000GP_hm3_Phase3_20k.bim" ] || [ ! -s "whole_genome_${ancestry}_1000GP_hm3_Phase3_20k.fam" ]; then
  ../../plink --bfile pre_merge_${ancestry}_1000GP_hm3_Phase3_chr1_20k.controls \
  --merge-list merge-list.txt \
  --make-bed \
  --allow-no-sex \
  --out whole_genome_${ancestry}_1000GP_hm3_Phase3_20k
fi

../../.././plink2 --bfile whole_genome_${ancestry}_1000GP_hm3_Phase3_20k --freq --out WG_af

echo "done with sim ${sim}"
