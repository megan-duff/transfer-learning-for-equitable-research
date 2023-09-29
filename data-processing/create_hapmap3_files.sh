#!/bin/bash

#pull second column (rsIDs) from .map file, remove any duplicate SNP IDs, and save under new file name for each ancestry group
awk '{print $2}' hapmap3_r1_b36_fwd.ASW.qc.poly.recode.map | sort -u > ASW_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.CEU.qc.poly.recode.map | sort -u > CEU_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map | sort -u > CHB_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.MEX.qc.poly.recode.map | sort -u > MEX_hapmap3_variants.txt
awk '{print $2}' hapmap3_r1_b36_fwd.TSI.qc.poly.recode.map | sort -u > TSI_hapmap3_variants.txt

#combine all hapmap3 snps files, pull second column (rsIDs) from .map file, remove any duplicate SNP IDs, and save under new file name 
awk '{print $2}' hapmap3_r1_b36_fwd.ASW.qc.poly.recode.map hapmap3_r1_b36_fwd.CEU.qc.poly.recode.map hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map hapmap3_r1_b36_fwd.CHB.qc.poly.recode.map hapmap3_r1_b36_fwd.MEX.qc.poly.recode.map hapmap3_r1_b36_fwd.TSI.qc.poly.recode.map | sort -u > merged_hapmap3_variants.txt
