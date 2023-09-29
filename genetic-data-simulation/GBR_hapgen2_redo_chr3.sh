#!/bin/bash
#SBATCH --job-name=GBR_HM3_GW
#SBATCH --output=GBR_HM3_GW_chr3
#SBATCH --error=GBR_HM3_GW_chr3_err
#SBATCH -n 16
#SBATCH -p math-alderaan
#SBATCH --array=232,245,275

#Breaks up 400 simulations of 22 chromosomes into 8800 jobs to be submitted at once

chr=3
sim_num=${SLURM_ARRAY_TASK_ID}

#Create dummy disease locus for HAPGEN2 because even though I am simulating under their control model it won't run wi$

dummyDL1=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr1_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL2=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr2_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL3=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr3_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL4=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr4_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL5=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr5_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL6=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr6_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL7=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr7_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL8=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr8_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL9=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr9_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL10=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr10_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL11=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr11_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL12=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr12_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL13=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr13_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL14=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr14_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL15=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr15_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL16=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr16_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL17=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr17_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL18=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr18_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL19=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr19_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL20=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr20_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL21=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr21_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL22=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr22_1000G_phase3_HM3.legend|head -40|tail -39)

dummy=dummyDL${chr}
dummy2=(`eval echo $dummy`)
dummy3=(`echo "${!dummy2}"`)

### simulate with HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html for details on HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html for details on effective sample size es$
### see https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html for details on files used for 1000G Phase 3 data on $

#Simulates 120000 GBR cases for each chromosome

./hapgen2 -m /storage/math/projects/duffme_gensim/Recom_Rates/GBR/GBR-${chr}-final.txt \
        -l /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr${chr}_1000G_phase3_HM3.legend \
        -h /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr${chr}_1000G_phase3_HM3.haps \
        -o /scratch/duffme_gensim/Simulations/GBR/sim_${sim_num}/GBR_1000GP_hm3_Phase3_chr${chr}_120k \
        -Ne 11418 \
        -n 120000 0 \
	-dl ${dummy3[0]} 1 1 1 ${dummy3[1]} 1 1 1 ${dummy3[2]} 1 1 1 ${dummy3[3]} 1 1 1 ${dummy3[4]} 1 1 1 ${dummy3[5]} 1 1 1 ${dummy3[6]} 1 1 1 ${dummy3[7]} 1 1 1 ${dummy3[8]} 1 1 1 ${dummy3[9]} 1 1 1 ${dummy3[10]} 1 1 1 ${dummy3[11]} 1 1 1 ${dummy3[12]} 1 1 1 ${dummy3[13]} 1 1 1 ${dummy3[14]} 1 1 1 ${dummy3[15]} 1 1 1 ${dummy3[16]} 1 1 1 ${dummy3[17]} 1 1 1 ${dummy3[18]} 1 1 1 ${dummy3[19]} 1 1 1 ${dummy3[20]} 1 1 1 ${dummy3[21]} 1 1 1 ${dummy3[22]} 1 1 1 ${dummy3[23]} 1 1 1 ${dummy3[24]} 1 1 1 ${dummy3[25]} 1 1 1 ${dummy3[26]} 1 1 1  ${dummy3[27]} 1 1 1 ${dummy3[28]} 1 1 1 ${dummy3[29]} 1 1 1 ${dummy3[30]} 1 1 1 ${dummy3[31]} 1 1 1 ${dummy3[32]} 1 1 1 ${dummy3[33]} 1 1 1 ${dummy3[34]} 1 1 1 ${dummy3[35]} 1 1 1 ${dummy3[36]} 1 1 1 ${dummy3[37]} 1 1 1 ${dummy3[38]} 1 1 1 ${dummy3[39]} 1 1 1 \
	-no_haps_output

if ! grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/GBR/sim_${sim_num}/GBR_1000GP_hm3_Phase3_chr${chr}_120k.summary
then
	echo chr: ${chr} sim: ${sim_num} >> GBR_failed_sims.txt
fi

if grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/GBR/sim_${sim_num}/GBR_1000GP_hm3_Phase3_chr${chr}_120k.summary
then
cd /scratch/duffme_gensim/Simulations/GBR/sim_${sim_num}
rm *cases*
for chr in {1..22};
do
../../.././plink --data GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls\
		--make-bed\
		--oxford-single-chr ${chr}\
		--out GBR_1000GP_hm3_Phase3_chr${chr}_120k.controls;
done
fi
