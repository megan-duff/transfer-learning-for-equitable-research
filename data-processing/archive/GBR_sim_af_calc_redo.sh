#!/bin/bash
#SBATCH --job-name=GBR_af_redo
#SBATCH --output=GBR_af_redo_%a.out
#SBATCH --error=GBR_af_redo_%a.err
#SBATCH -n 10
#SBATCH -p math-alderaan
#SBATCH --array=1-400

sim=${SLURM_ARRAY_TASK_ID}

cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim}

if [ ! -s "whole_genome_GBR_1000GP_hm3_Phase3_20k.bed" ] || [ ! -s "whole_genome_GBR_1000GP_hm3_Phase3_20k.bim" ] || [ ! -s "whole_genome_GBR_1000GP_hm3_Phase3_20k.fam" ]; then
    ../../plink --bfile subset_1_whole_genome_GBR_1000GP_hm3_Phase3_20k \
    --merge-list merge-list.txt \
    --make-bed \
    --allow-no-sex \
    --out whole_genome_GBR_1000GP_hm3_Phase3_20k
fi

if [ ! -s "WG_af.afreq" ]; then
    ../../.././plink2 --bfile whole_genome_GBR_1000GP_hm3_Phase3_20k --freq --out WG_af
fi

echo "done with sim ${sim}"
