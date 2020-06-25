#!/bin/bash
for i in /scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr*.traw; do
  Rscript SKAT_multithreaded_wrapper.R \
  --phenodata_path $1 \
  --covariates /scratch2/NSF_GWAS/phenodata/Covariates_DiamFinal_PhaseFinal_PCs.txt \
  --raw_file_path $i \
  --window_size 3000 \
  --window_shift 1000 \
  --output_dir /scratch2/NSF_GWAS/Results/SKAT/
done

