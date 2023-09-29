#!/bin/bash
#SBATCH --job-name=GBR_gw_files_redo
#SBATCH --output=GBR_gw_files_redo_%a.out
#SBATCH --output=GBR_gw_files_redo_%a.err
#SBATCH -n 15
#SBATCH -p math-alderaan
#SBATCH --array=347,372-375

sim=${SLURM_ARRAY_TASK_ID}
cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

rm merge-list.txt

for i in {2..22}; do
    echo "GBR_1000GP_hm3_Phase3_chr${i}_120k.controls" >> merge-list.txt
done

subset_1_file="subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k-merge.missnp"

if [ ! -s "$subset_1_file" ]; then
  ../../plink --bfile GBR_1000GP_hm3_Phase3_chr1_120k.controls \
  --merge-list merge-list.txt \
  --keep subset_1.txt \
  --make-bed \
  --out subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k
fi

for chr in {1..22}; do 

  bed_file="GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.bed"
  bim_file="GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.bim"
  fam_file="GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp.fam"
  
  if [ ! -s "$bed_file" ] || [ ! -s "$bim_file" ] || [ ! -s "$fam_file" ]; then
  
    ../../plink --bfile GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls\
    --make-bed\
    --exclude subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k-merge.missnp\
    --out GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp;
  fi;
  
done 


rm merge-list.txt

for chr in {2..22}; do
    echo "GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls_no_missnp" >> merge-list.txt
done

#If file does not exist in current directory, then re-run code to create the file 

for subset in {1..6}; do

  bed_file="subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k.bed"
  bim_file="subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k.bim"
  fam_file="subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k.fam"

  if [ ! -s "$bed_file" ] || [ ! -s "$bim_file" ] || [ ! -s "$fam_file" ]; then
    echo "Sim ${sim}: File doesn't exist for subset ${subset}"
    ../../plink --bfile GBR_1000GP_hm3_Phase3_chr1_120k.controls_no_missnp\
    --merge-list merge-list.txt \
    --keep subset_${subset}.txt \
    --make-bed \
    --allow-no-sex\
    --out subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k
    fi;

done

echo "Done with sim ${sim}"
