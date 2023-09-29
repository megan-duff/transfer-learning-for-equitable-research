#!/bin/bash
#SBATCH --job-name=CEU_HM3_GW_redo
#SBATCH --output=CEU_HM3_GW_redo
#SBATCH --error=CEU_HM3_GW_redo_err
#SBATCH -n 16
#SBATCH -p math-alderaan
#SBATCH --array=7504,7566,7700,7722,7727,7753,7754,7755,7756,7757,7758,7759,7760,7761,7762,7763,7764,7765,7767,7768,7769,7770,7771,7772,7773,7774,7776,7777,7778,7779,7780,7781,7782,7783,7784,7785,7786,7787,7788,7789,7790,7791,7792,7793,7794,7795,7796,7797,7798,7799,7800,7801,7802,7803,7804,7805,7806,7807,7808,7809,7810,7811,7812,7813,7814,7815,7816,7817,7818,7819,7820,7821,7822,7823,7824,7825,7826,7827,7828,7829,7830,7831,7832,7833,7834,7835,7836,7837,7838,7839,7840,7841,7842,7843,7844,7845,7846,7847,7848,7849,7850,7851,7852,7853,7854,7855,7856,7857,7858,7859,7860,7861,7862,7863,7864,7865,7866,7867,7868,7869,7871,7872,7873,7874,7875,7876,7877,7878,7879,7880,7881,7882,7883,7884,7885,7886,7887,7888,7889,7890,7891,7892,7893,7894,7895,7896,7897,7898,7899,7900,7901,7902,7903,7904,7905,7906,7907,7908,7909,7910,7911,7912,7913,7914,7915,7916,7917,7918,7919,7920,7921,7922,7925,7926,7927,7929,7930,7931,7932,7933,7934,7935,7936,7937,7938,7939,7940,7941,7942,7943,7944,7945,7946,7947,7948,7949,7950,7951,7952,7953,7954,7955,7956,7957,7958,7959,7960,7961,7962,7963,7964,7965,7966,7967,7968,7969,7970,7971,7972,7973,7974,7975,7976,7977,7978,7979,7980,7981,7982,7984,7985,7986,7987,7988,7989,7990,7991,7992,7993,7994,7995,7996,7997,7998,7999,8000,8001,8002,8003,8004,8005,8006,8007,8008,8009,8010,8011,8012,8013,8014,8015,8016,8017,8018,8019,8020,8021,8022,8023,8024,8025,8026,8029,8030,8031,8032,8033,8034,8035,8036,8037,8038,8039,8040,8041,8042,8043,8044,8045,8046,8047,8048,8049,8050,8051,8052,8053,8054,8055,8056,8057,8058,8059,8060,8061,8062,8063,8064,8065,8066,8067,8068,8069,8070,8071,8072,8073,8074,8075,8076,8077,8078,8079,8080

#Breaks up 400 simulations of 22 chromosomes into 8400 jobs to be submitted at once

chr=$(expr ${SLURM_ARRAY_TASK_ID} % 22 + 1)
((sim_num=${SLURM_ARRAY_TASK_ID}/22 + 1))

#Create dummy disease locus for HAPGEN2 because even though I am simulating under their control model it won't run wi$

dummyDL1=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr1_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL2=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr2_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL3=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr3_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL4=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr4_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL5=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr5_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL6=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr6_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL7=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr7_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL8=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr8_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL9=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr9_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL10=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr10_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL11=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr11_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL12=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr12_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL13=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr13_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL14=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr14_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL15=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr15_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL16=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr16_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL17=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr17_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL18=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr18_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL19=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr19_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL20=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr20_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL21=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr21_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL22=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr22_1000G_phase3_HM3.legend|head -40|tail -39)

dummy=dummyDL${chr}
dummy2=(`eval echo $dummy`)
dummy3=(`echo "${!dummy2}"`)

### simulate with HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html for details on HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html for details on effective sample size es$
### see https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html for details on files used for 1000G Phase 3 data on $

#Simulates 20000 CEU cases for each chromosome

./hapgen2 -m /storage/math/projects/duffme_gensim/Recom_Rates/CEU/CEU-${chr}-final.txt \
        -l /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr${chr}_1000G_phase3_HM3.legend \
        -h /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr${chr}_1000G_phase3_HM3.haps \
        -o /scratch/duffme_gensim/Simulations/CEU/sim_${sim_num}/CEU_1000GP_hm3_Phase3_chr${chr}_20k \
        -Ne 11418 \
        -n 20000 0 \
	-dl ${dummy3[0]} 1 1 1 ${dummy3[1]} 1 1 1 ${dummy3[2]} 1 1 1 ${dummy3[3]} 1 1 1 ${dummy3[4]} 1 1 1 ${dummy3[5]} 1 1 1 ${dummy3[6]} 1 1 1 ${dummy3[7]} 1 1 1 ${dummy3[8]} 1 1 1 ${dummy3[9]} 1 1 1 ${dummy3[10]} 1 1 1 ${dummy3[11]} 1 1 1 ${dummy3[12]} 1 1 1 ${dummy3[13]} 1 1 1 ${dummy3[14]} 1 1 1 ${dummy3[15]} 1 1 1 ${dummy3[16]} 1 1 1 ${dummy3[17]} 1 1 1 ${dummy3[18]} 1 1 1 ${dummy3[19]} 1 1 1 ${dummy3[20]} 1 1 1 ${dummy3[21]} 1 1 1 ${dummy3[22]} 1 1 1 ${dummy3[23]} 1 1 1 ${dummy3[24]} 1 1 1 ${dummy3[25]} 1 1 1 ${dummy3[26]} 1 1 1  ${dummy3[27]} 1 1 1 ${dummy3[28]} 1 1 1 ${dummy3[29]} 1 1 1 ${dummy3[30]} 1 1 1 ${dummy3[31]} 1 1 1 ${dummy3[32]} 1 1 1 ${dummy3[33]} 1 1 1 ${dummy3[34]} 1 1 1 ${dummy3[35]} 1 1 1 ${dummy3[36]} 1 1 1 ${dummy3[37]} 1 1 1 ${dummy3[38]} 1 1 1 ${dummy3[39]} 1 1 1 \
	-no_haps_output

if ! grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/CEU/sim_${sim_num}/CEU_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
	echo chr: ${chr} sim: ${sim_num} >> CEU_failed_sims.txt
fi

if grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/CEU/sim_${sim_num}/CEU_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
cd /scratch/duffme_gensim/Simulations/CEU/sim_${sim_num}
rm *cases*
for chr in {1..22};
do
../../.././plink --data CEU_1000GP_hm3_Phase3_chr${chr}_20k.controls\
		--make-bed\
		--allow-no-sex\
		--oxford-single-chr ${chr}\
		--out CEU_1000GP_hm3_Phase3_chr${chr}_20k.controls;
done
fi
