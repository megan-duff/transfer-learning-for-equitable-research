## Data Processing and Genetic Data Simulation README

This read me outlines to steps taken to download and process the reference data and then simulate genetic data. This simulation study will be used to compare various transfer learning methods against baseline methods (methods previously used for genomic prediction). I will simulate 120,000 individuals of GBR (British in England and Scotland) genetic ancestry and 20,000 individuals of each of the following ancestries:YRI (Yorban in Africa), CHB (Chinese in Bejing, China), TSI (Toscani in Italy), and CEU (Northern and Western European in Utah). 

These individuals will be simulated using the software HAPGEN2 (https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html), where the reference population is from the corresponding genetic ancestry from the 1000 Genomes Project (http://www.internationalgenome.org/). For this reserch, Build37 will be used which has a total of 77,818,332 biallelic SNPs. From this set of SNPs, only the variants identified in the HapMap3 project will be simulated. 

Overall, this chapter outlines the bioinformatics pipeline to simulate the appropriate data used for method comparison. This read me / code is associated with Chapter 3 of my dissertation. 


### Prepare Reference Data

To prepare the reference data, the first step is to download the appropriate data from the following sites. 

*Genetic map (recombination rates) for each subpopulation: ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates

*1000G Build 37 Data: https://www.cog-genomics.org/plink/2.0/resources

*HapMap3 data: https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/

To download 1000G Build 37 Data, both whole-genome and chromosome specific files are downloaded. Chromosome specific files will be used for HAPGEN2 where whole-genome files will be used for analyses, such as principal components analysis. These files are provided on https://www.cog-genomics.org/plink/2.0/resources.
```{bash, eval=F}
#Download whole-genome 1000G project data
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1
wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1
wget https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1

#Download HapMap3 Data
wget "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2"

#Download recom rate files and untar files
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CHB_omni_recombination_20130507.tar
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CEU_omni_recombination_20130507.tar
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/GBR_omni_recombination_20130507.tar
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CLM_omni_recombination_20130507.tar
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/TSI_omni_recombination_20130507.tar
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/YRI_omni_recombination_20130507.tar

tar -xf CHB_omni_recombination_20130507.tar
tar -xf CEU_omni_recombination_20130507.tar
tar -xf GBR_omni_recombination_20130507.tar
tar -xf CLM_omni_recombination_20130507.tar
tar -xf TSI_omni_recombination_20130507.tar
tar -xf YRI_omni_recombination_20130507.tar

#Download chromosome specific 1000G project data
wget https://www.dropbox.com/s/5gbbt3z42z652xt/chr1_phase3.pgen.zst?dl=1 &
wget https://www.dropbox.com/s/6qlhq2mdawa27f2/chr2_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/nimp2m8z4iqhqh8/chr3_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/zubmyfwsyhnbawb/chr4_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/fa4r8skrw2mqfvo/chr5_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/la56okhu4p2ul18/chr6_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/gknbjkflrngxt3n/chr7_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/zshg2k7jod3w5e1/chr8_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/8nrprvfamr5fwf8/chr9_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/6tm2rf6zfgelwod/chr10_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/5syf1fjq2u95u8r/chr11_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/hg36svzdoynecfm/chr12_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/fzx42kpxeg4ktkg/chr13_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/xby4sgwwwt333bl/chr14_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/qloe20hm8t0hfqg/chr15_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/j6cuuc8q1jzwj9y/chr16_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/gga4fevcvlvb92k/chr17_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/n4qbg7j69ldtex3/chr18_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/ib5wiwz71yk1qif/chr19_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/f9nqaps0ddsguxe/chr20_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/kh120unohlk2e37/chr21_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/w9wwua4pe9em280/chr22_phase3.pgen.zst?dl=1 & 
wget https://www.dropbox.com/s/hw5qsh45rofqkxg/chr1_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/bof8v3odxtd8ihm/chr2_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/cjj2c6kafulzg4e/chr3_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/vmic0acyuru2ojl/chr4_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/dfgxjbl5j0dlony/chr5_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/iko2dn5565hvyqn/chr6_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/m5s6kwbayoi5266/chr7_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/rrs3zx5spkvuyjo/chr8_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/wba357vnqjibkmd/chr9_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/8mu9qf45wcmattv/chr10_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/05jypoy5mvlb2va/chr11_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/s26n11e6yq510hc/chr12_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/cbxbkczjb2s69rz/chr13_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/vm7dqyljembn62g/chr14_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/5zxfjqfyadzbidg/chr15_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/5ynstl55kcahpmr/chr16_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/ufruz7gvrjx7f0e/chr17_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/1oqweb2dhqcwuo4/chr18_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/0x8bglz4hkeoaxy/chr19_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/mk7clb8pfyghikl/chr20_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/uydnboqp4rlfryf/chr21_phase3.pvar.zst?dl=1 & 
wget https://www.dropbox.com/s/c4kqgoc93sir2g5/chr22_phase3.pvar.zst?dl=1 & 

#clean up file format / file names
.././plink2 --zst-decompress 'all_phase3.pgen.zst?dl=1' > all_phase3.pgen
.././plink2 --zst-decompress 'all_phase3.pvar.zst?dl=1' > all_phase3.pvar

for chr in {1..22}; do .././plink2 --zst-decompress chr${chr}_phase3.pvar.zst?dl=1 > chr_${chr}.pvar; done

for chr in {1..22}; do .././plink2 --zst-decompress chr${chr}_phase3.pgen.zst?dl=1 > chr_${chr}.pgen; done
```

To create the appropriate HapMap3 SNP files, the following files are created:

* ASW list of HapMap3 variants : ASW_hapmap3_variants.txt
  * 1,536,247 variants 
* CEU list of HapMap3 variants : CEU_hapmap3_variants.txt
  * 1,403,896 variants
* CHB list of HapMap3 variants : CHB_hapmap3_variants.txt
  * 1,311,113 variants 
* MEX list of HapMap3 variants : MEX_hapmap3_variants.txt
  * 1,430,334 variants
* TSI list of HapMap3 variants : TSI_hapmap3_variants.txt
  * 1,393,925 variants 
* ASW, CEU, CHB, MEX, and TSI HapMap3 variants : merged_hapmap3_variants.txt
  * 1,603,825 variants 
  
```{bash, eval =F}
#!/bin/bash
#pull second column (rsIDs) from .map file, remove any duplicate SNP IDs, and save under new file name for each ancestry group
awk '{print $2}' hapmap3_r1_b36_fwd.ASW.qc.poly.recode.map | sort -u > ASW_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.CEU.qc.poly.recode.map | sort -u > CEU_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map | sort -u > CHB_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.MEX.qc.poly.recode.map | sort -u > MEX_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.TSI.qc.poly.recode.map | sort -u > TSI_hapmap3_variants.txt

#combine all hapmap3 snps files, pull second column (rsIDs) from .map file, remove any duplicate SNP IDs, and save under new file name 
awk '{print $2}' hapmap3_r1_b36_fwd.ASW.qc.poly.recode.map hapmap3_r1_b36_fwd.CEU.qc.poly.recode.map hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map hapmap3_r1_b36_fwd.MEX.qc.poly.recode.map hapmap3_r1_b36_fwd.TSI.qc.poly.recode.map | sort -u > merged_hapmap3_variants.txt

```

Note: Originally the plan was to simulate the five target populations: ASW, CEU, CHB, MEX, and TSI. Due to HAPGEN2's inability to accurately simulate admixed populations, we cannot simulate ASW and MEX. Therefore, YRI was added to the list of target populations that will be simulated instead. We are still using the HapMap3 variants for the original five populations. 

#### Data Cleaning for WG files and for Chromosome specific files

Next, the 1000G Project data was quality controlled. Only bi-allelic HAPMAP 3 SNPs were kept. Note, when HapMap3 SNPs are referenced, this is referring to the combined set of HapMap3 SNPs across the 5 population groups (file: merged_hapmap3_variants.txt). Additionally, first and second degree relatives were removed from the set. Note, plink2 resources provide a list of samples that have a first or second degree relative (file: 'deg2_phase3.king.cutoff.out.id?dl=1'), which is used in plink to remove said individuals. 

```{bash, eval = F}
#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out WG_1000G_phase3

#Keep only autosomal SNPs, extract hapmap3 variants, remove first and second degree relatives, and maf filter of 0.05 for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --maf 0.05 --out WG_1000G_phase3_common

#Keep only autosomal SNPs, extract hapmap3 variants, remove first and second degree relatives, maf filter of 0.05, and LD prune for whole genome files
../.././plink2 --make-bed --pfile all_phase3 --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --indep-pairwise 50 5 0.5 --maf 0.05 --out WG_1000G_phase3_HM3_common_ldpruned

#rename .psam file to have the same prefix as the .pgen files (this makes it easier to input into plink)
for chr in {1..22}; do cp all_phase3.psam chr_${chr}.psam; done

#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for chromosome specific files
for chr in {1..22}; do 
../.././plink2 --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do mv chr_${chr}.psam chr${chr}_1000G_phase3.psam ; done

#Keep only autosomal SNPs, extract hapmap3 variants, and remove first and second degree relatives for chromosome specific files & save the file as .haps/.legend for HAPGEN2
for chr in {1..22}; do 
../.././plink2 --export hapslegend --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --out chr${chr}_1000G_phase3_HM3; done
```

After removing related individuals and subseting to the merged HapMap3 variants, the files contained 2,490 samples and 1,520,315 genome wide variants. 

#### Create Reference Files for HAPGEN2

Each population to be simulated (GBR, YRI, CEU, TSI, CHB) will need a corresponding reference file. 

First, create text files with the corresponding IDs for each population group. The all_phase3.psam file contains the population and super population group for each sample in the 1000G project data. This will be used to create files that contain the IDs that correspond to each of the six ancestry groups. These files will then be loaded into plink to subset the full 1000G project data into ancestry-specific files for each of the six ancestry groups. 
To create the text files of IDs: 
```{R, eval = F}
singularity shell /storage/singularity/mixtures.sif
R
#Read in data that contains 
data = read.csv('all_phase3.psam', header = T, sep="\t")

GBR = subset(data, data$Population == "GBR")
GBR_ids = GBR[,1]
write.table(GBR_ids, "GBR_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

YRI = subset(data, data$Population == "YRI")
YRI_ids = YRI[,1]
write.table(YRI_ids, "YRI_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

CEU = subset(data, data$Population == "CEU")
CEU_ids = CEU[,1]
write.table(CEU_ids, "CEU_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

CHB = subset(data, data$Population == "CHB")
CHB_ids = CHB[,1]
write.table(CHB_ids, "CHB_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

TSI = subset(data, data$Population == "TSI")
TSI_ids = TSI[,1]
write.table(TSI_ids, "TSI_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

CLM = subset(data, data$Population == "CLM")
CLM_ids = CLM[,1]
write.table(CLM_ids, "CLM_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

Next, subset 1000G_HM3 files to only include samples corresponding to each ancestry group:
```{bash, eval = F}
#Create .haps/.legend for HAPGEN2 for each ancestry group
# Note, the 1000G original files are fully phased and so the command -export hapslegend works 
# If you try to use the subsetted file that was used for analysis earlier, it will not work 

for chr in {1..22}; do 
../.././plink2 --export hapslegend --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep GBR_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --export hapslegend --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep CEU_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --export hapslegend --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep CHB_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --export hapslegend --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep CLM_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CLM/CLM_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --export hapslegend --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep TSI_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --export hapslegend --make-bed --pfile chr_${chr} --max-alleles 2 --autosome --memory 8000 --extract /storage/math/projects/duffme_gensim/1000G_Data/HapMap3/merged_hapmap3_variants.txt --snps-only --remove 'deg2_phase3.king.cutoff.out.id?dl=1' --keep YRI_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/YRI/YRI_chr${chr}_1000G_phase3_HM3; done

#Create plink format files for analysis of each ancestry group
for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep GBR_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/GBR/GBR_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep ASW_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/ASW/ASW_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep CEU_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CEU/CEU_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep CHB_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/CHB/CHB_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep MXL_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/MXL/MXL_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep TSI_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr${chr}_1000G_phase3_HM3; done

for chr in {1..22}; do 
../.././plink2 --make-pgen --bfile chr${chr}_1000G_phase3_HM3 --keep YRI_ids.txt --out /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/YRI/YRI_chr${chr}_1000G_phase3_HM3; done

```

### Genetic Data Simulation: HAPGEN2 Simulation Scripts
Now that the reference files have been created and quality controlled, HAPGEN2 is used to simulate genetic data through the following bash scripts.

The inputs for HAPGEN2 are as follows:

* -m: the genetic map or recombination rate file (ancestry specific)
* -l: legend file for samples.
* -h: haplotype file for samples.
* -o: output file to be created. 
* -Ne: effective sample size (ancestry specific). Using estimates from https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html. 
* -n: Number of desired simuulated control and case individuals. All individuals will be simulated as controls in this workflow. 
* dl: disease or causal loci. Note, since we are simulating controls this does not matter but HAPGEN2 will not run without this parameter. Therefore, dummy causal loci are selected from each chromosome as a placeholder. 

The following script to simulate the GBR individuals is given as an example. Note, the other ancestry group simulation scripts follow the same template but vary in  their corresponding ancestry input files, effective sample size estimate, and number of simulated individuals. Recall, the simulation framework simulates 120,000 GBR individuals and 20,000 ASW, TSI, CEU, CHB, and MXL individuals. Additionally, the following estimates of effective sample size are used: 

* GBR: 11,418
* TSI: 11,418
* CEU: 11,418
* YRI: 17,469
* CHB: 14,269

GBR script: 
```{bash GBR script, eval = F}
#!/bin/bash
#SBATCH --job-name=GBR_HM3_GW
#SBATCH --output=GBR_HM3_GW
#SBATCH --error=GBR_HM3_GW_err
#SBATCH -n 16
#SBATCH -p math-alderaan
#SBATCH --array 0-8799

#Breaks up 400 simulations of 22 chromosomes into 8800 jobs to be submitted at once

chr=$(expr ${SLURM_ARRAY_TASK_ID} % 22 + 1)
((sim_num=${SLURM_ARRAY_TASK_ID}/22+1))

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
```

Similar scripts were used to simulate TSI, CEU, CHB, and YRI descent individuals. 
