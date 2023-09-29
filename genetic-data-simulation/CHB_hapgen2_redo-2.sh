#!/bin/bash
#SBATCH --job-name=CHB_HM3_GW
#SBATCH --output=CHB_HM3_GW
#SBATCH --error=CHB_HM3_GW_err
#SBATCH -n 16
#SBATCH -p math-alderaan
#SBATCH --array 23,155,419,506,595,705,1210,1211,1277,1347,1364,1474,1848,1915,1958,1959,2047,2077,2202,2487,2577,2596,2927,3344,3388,3389,3411,3432,3433,3498,3630,3851,4115,4136,4270,4423,4466,4621,4642,4643,5390,5478,5567,5677,5744,6138,6227,6381,6535,6578,7018,7128,7129,7172,7196,7569,7590,7634,7656,7723,7745,7810,7855,7965,8031,8382,8383,8470,8471,8514,8626,8627,8630,8631

#Breaks up 400 simulations of 22 chromosomes into 8400 jobs to be submitted at once

chr=$(expr ${SLURM_ARRAY_TASK_ID} % 22 + 1)
((sim_num=${SLURM_ARRAY_TASK_ID}/22 + 1))

#Create dummy disease locus for HAPGEN2 because even though I am simulating under their control model it won't run wi$

dummyDL1=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr1_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL2=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr2_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL3=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr3_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL4=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr4_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL5=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr5_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL6=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr6_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL7=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr7_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL8=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr8_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL9=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr9_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL10=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr10_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL11=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr11_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL12=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr12_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL13=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr13_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL14=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr14_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL15=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr15_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL16=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr16_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL17=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr17_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL18=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr18_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL19=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr19_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL20=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr20_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL21=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr21_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL22=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr22_1000G_phase3_HM3.legend|head -40|tail -39)

dummy=dummyDL${chr}
dummy2=(`eval echo $dummy`)
dummy3=(`echo "${!dummy2}"`)

### simulate with HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html for details on HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html for details on effective sample size 
### see https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html for details on files used for 1000G Phase 3 data 

#Simulates 20000 CHB cases for each chromosome

./hapgen2 -m /storage/math/projects/duffme_gensim/Recom_Rates/CHB/CHB-${chr}-final.txt \
        -l /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr${chr}_1000G_phase3_HM3.legend \
        -h /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr${chr}_1000G_phase3_HM3.haps \
        -o /scratch/duffme_gensim/Simulations/CHB/sim_${sim_num}/CHB_1000GP_hm3_Phase3_chr${chr}_20k \
        -Ne 14269 \
        -n 20000 0 \
	-dl ${dummy3[0]} 1 1 1 ${dummy3[1]} 1 1 1 ${dummy3[2]} 1 1 1 ${dummy3[3]} 1 1 1 ${dummy3[4]} 1 1 1 ${dummy3[5]} 1 1 1 ${dummy3[6]} 1 1 1 ${dummy3[7]} 1 1 1 ${dummy3[8]} 1 1 1 ${dummy3[9]} 1 1 1 ${dummy3[10]} 1 1 1 ${dummy3[11]} 1 1 1 ${dummy3[12]} 1 1 1 ${dummy3[13]} 1 1 1 ${dummy3[14]} 1 1 1 ${dummy3[15]} 1 1 1 ${dummy3[16]} 1 1 1 ${dummy3[17]} 1 1 1 ${dummy3[18]} 1 1 1 ${dummy3[19]} 1 1 1 ${dummy3[20]} 1 1 1 ${dummy3[21]} 1 1 1 ${dummy3[22]} 1 1 1 ${dummy3[23]} 1 1 1 ${dummy3[24]} 1 1 1 ${dummy3[25]} 1 1 1 ${dummy3[26]} 1 1 1  ${dummy3[27]} 1 1 1 ${dummy3[28]} 1 1 1 ${dummy3[29]} 1 1 1 ${dummy3[30]} 1 1 1 ${dummy3[31]} 1 1 1 ${dummy3[32]} 1 1 1 ${dummy3[33]} 1 1 1 ${dummy3[34]} 1 1 1 ${dummy3[35]} 1 1 1 ${dummy3[36]} 1 1 1 ${dummy3[37]} 1 1 1 ${dummy3[38]} 1 1 1 ${dummy3[39]} 1 1 1 \
	-no_haps_output

if ! grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/CHB/sim_${sim_num}/CHB_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
	echo chr: ${chr} sim: ${sim_num} >> CHB_failed_sims.txt
fi

if grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/CHB/sim_${sim_num}/CHB_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
cd /scratch/duffme_gensim/Simulations/CHB/sim_${sim_num}
rm *cases*
for chr in {1..22};
do
../../.././plink --data CHB_1000GP_hm3_Phase3_chr${chr}_20k.controls\
		--make-bed\
		--allow-no-sex\
		--oxford-single-chr ${chr}\
		--out CHB_1000GP_hm3_Phase3_chr${chr}_20k.controls;
done
fi
