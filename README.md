# merge_multistudy_data
This repo is to merge genotype (DNA) and phenotype data (gene expression feature counts) of multiple datasets having common cell_types

Steps:
1. renames feature counts files to be unique
E.g: feature_counts_merged.txt -> CAP_feature_counts_merged.txt

2. Index/Merge/Filter the provided VCf files

3. selects the QC passed samples and merges feature_counts files
  - Merges the samples metadata files
  - Selects the samples which passed the QC and have the provided "common_cell_type" with "naive" condition
  - selects the genotype_ids to be extracted from merged_vcf file and writes them into the file
  - merges the Feature counts files (no sample selection here) and checks if all of them has the same number of genes and if same sample_id occurs in different studies
  
 4. extracts the genotype samples (which we selected in step 3) from the merged VCF file
