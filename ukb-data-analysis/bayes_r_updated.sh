
# Download GCTB
wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.1_Linux.zip
unzip gctb_2.5.1_Linux.zip

# Download Plink 
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip plink_linux_x86_64_20231211.zip

# Set variables
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

dx download -f Thesis:/Genotype_Data/target_clump_snps/${ancestry}_whole_genome_maf_unrelated_pheno_${pheno_code}_1neg2.clumps

awk '{print $3}' ${ancestry}_whole_genome_maf_unrelated_pheno_${pheno_code}_1neg2.clumps > pheno_${pheno_code}_rs_ids.txt

# Run Bayes R via GCTB 
./gctb_2.5.1_Linux/gctb --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train --extract pheno_${pheno_code}_rs_ids.txt --pheno ${ancestry}_corrected_predicted_${pheno_name}_train_formatted.txt --bayes R --chain-length 12000 --burn-in 2000 --out ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2

# Compute train set bayes R score
./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_train --score ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.snpRes 2 5 8 header sum center --out train_pred_${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2

# Compute validation set bayes R score
./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_val --score ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.snpRes 2 5 8 header sum center --out val_pred_${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2

# Compute test set bayes R score
./plink --bfile ${ancestry}_whole_genome_maf_unrelated_clump_snps_1e-2_test --score ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.snpRes 2 5 8 header sum center --out test_pred_${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2

mkdir ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
mv ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.covRes ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
mv ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.snpRes ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
mv ${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.parRes ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
mv *.log ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
mv *_${ancestry}_pheno_${pheno_code}_bayes_r_snps_clump_snps_1e-2.profile ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2
dx upload -r ${ancestry}_pheno_${pheno_code}_bayes_r_clump_snps_1e-2