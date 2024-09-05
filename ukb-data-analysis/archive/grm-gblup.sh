dx download Thesis:/Software/gcta-1.94.1-linux-kernel-3-x86_64.zip

unzip gcta-1.94.1-linux-kernel-3-x86_64.zip

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip

unzip plink_linux_x86_64_20231211.zip

ancestry=$1
pheno_name=$2
pheno_code=$3

# Download train/val/test plink files
for set in 'train' 'val' 'test'; do 
dx download -f Thesis:/Genotype_Data/${ancestry}/target_model_files/${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_${set}.bed
dx download -f Thesis:/Genotype_Data/${ancestry}/target_model_files/${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_${set}.bim
dx download -f Thesis:/Genotype_Data/${ancestry}/target_model_files/${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_${set}.fam;
done

# Download phenotype files 
dx download -f Thesis:/Data/phenotypes/${ancestry}_corrected_predicted_${pheno_name}_train.txt
dx download -f Thesis:/Data/phenotypes/${ancestry}_corrected_predicted_${pheno_name}_val.txt
dx download -f Thesis:/Data/phenotypes/${ancestry}_corrected_predicted_${pheno_name}_test.txt

# Format phenotype files to succesfully match plink format
# Plink files are in format: FID, IID where FID and IID are the same
# Phenotype files are in format: IID Phenotype
# Therefore need to duplicate IID column in phenotype file
input_file="${ancestry}_corrected_predicted_${pheno_name}_train.txt"
output_file="${ancestry}_corrected_predicted_${pheno_name}_train_formatted.txt"
awk 'NR>1 {print $1, $1, $2}' $input_file > $output_file

dx download -f Thesis${ancestry}_whole_genome_maf_unrelated_pheno_${pheno_code}_1neg2.clumps

awk '{print $3}' ${ancestry}_whole_genome_maf_unrelated_pheno_${pheno_code}_1neg2.clumps > pheno_${pheno_code}_rs_ids.txt

file_name="${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2" &&

for i in {1..10};
do ./gcta-1.94.1-linux-kernel-3-x86_64/gcta64 --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train --extract pheno_${pheno_code}_rs_ids.txt --make-grm-part 10 ${i} --out ${file_name};
done &&

cat ${file_name}.part_10_*.grm.id > ${file_name}.grm.id &&
cat ${file_name}.part_10_*.grm.bin > ${file_name}.grm.bin &&
cat ${file_name}.part_10_*.grm.N.bin > ${file_name}.grm.N.bin &&

./gcta-1.94.1-linux-kernel-3-x86_64/gcta64 --reml --grm ${file_name} --pheno ${ancestry}_corrected_predicted_${pheno_name}_train_formatted.txt --thread-num 64 --reml-pred-rand --out ${file_name}

./gcta-1.94.1-linux-kernel-3-x86_64/gcta64 --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train --extract pheno_${pheno_code}_rs_ids.txt --blup-snp ${file_name}.indi.blp --out ${file_name} &&

./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train --score ${file_name}.snp.blp sum 1 2 3 --out train_pred_${file_name} &&

./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_val --score ${file_name}.snp.blp sum 1 2 3 --out val_pred_${file_name} &&

./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_test --score ${file_name}.snp.blp sum 1 2 3 --out test_pred_${file_name} && 

mkdir ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2 &&
mv ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2.indi.blp ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2 &&
mv ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2.snp.blp ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2 &&
mv *_${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2.profile ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2 &&
dx upload -r ${ancestry}_pheno_${pheno_code}_gblup_clump_snps_1e-2
