#!/bin/bash
#SBATCH --job-name=GBR_gw_files
#SBATCH --output=GBR_gw_files_%j
#SBATCH --error=GBR_gw_files_err_%j
#SBATCH -n 10
#SBATCH -p math-alderaan
#SBATCH --array=11-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

rm subset_*.txt

### Break up 120k individuals into 6 sets of 20k people
for person in {0..19999}; do echo id1_${person} id2_${person} >> subset_1.txt; done
for person in {20000..39999}; do echo id1_${person} id2_${person} >> subset_2.txt; done
for person in {40000..59999}; do echo id1_${person} id2_${person} >> subset_3.txt; done
for person in {60000..79999}; do echo id1_${person} id2_${person} >> subset_4.txt; done
for person in {80000..99999}; do echo id1_${person} id2_${person} >> subset_5.txt; done
for person in {100000..119999}; do echo id1_${person} id2_${person} >> subset_6.txt; done

### Subset simulated genotype files to only include  specified samples

rm merge-list.txt

for i in {2..22}; do
    echo "GBR_1000GP_hm3_Phase3_chr${i}_120k.controls" >> merge-list.txt
done

echo "subset and merge"

../../plink --bfile GBR_1000GP_hm3_Phase3_chr1_120k.controls \
--merge-list merge-list.txt \
--keep subset_1.txt \
--make-bed \
--out subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k

echo "get rid of duplicate"

for chr in {1..22}; do
    ../../plink --bfile GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls\
    --make-bed\
    --exclude subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k-merge.missnp\
    --out GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp;
done

rm merge-list.txt

for chr in {2..22}; do
    echo "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp" >> merge-list.txt
done

echo "now real merge"

for subset in {1..6}; do
    ../../plink --bfile GBR_1000GP_hm3_Phase3_chr1_120k.controls_no_missnp \
    --merge-list merge-list.txt \
    --keep subset_${subset}.txt \
    --make-bed \
    --allow-no-sex\
    --out subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k;
done

for subset in {1..6}; do wc -l subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k.bim; done
