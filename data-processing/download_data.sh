#!/bin/bash

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
