#!/bin/bash
#SBATCH --job-name=GBR_af_redo
#SBATCH --output=GBR_af_redo_%a.out
#SBATCH --error=GBR_af_redo_%a.err
#SBATCH -n 10
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

rm merge-list.txt

for subset in {2..6}; do
    echo "subset_${subset}_whole_genome_GBR_1000GP_hm3_Phase3_20k" >> merge-list.txt
done

../../plink --bfile subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k\
    --merge-list merge-list.txt \
    --make-bed \
    --allow-no-sex\
    --out whole_genome_GBR_1000GP_hm3_Phase3_20k

../../.././plink2 --bfile whole_genome_GBR_1000GP_hm3_Phase3_20k --freq --out WG_af

echo "done with sim ${sim}"
