phenodata <- system.file("extdata",
                         "TDZ_shoot_area.plink.pheno",
                         package = "SKATMCMT")

# phenodata <- "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/shoot_5w.header.pheno"
#
# covariates <- "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt"

covariates <- system.file("extdata",
                          "poplar_PCs_covariates.tbt",
                          package = "SKATMCMT")

covariates <- system.file("extdata",
                          "poplar_PCs_covariates.tbt",
                          package = "SKATMCMT")

raw_file_path <- system.file("extdata",
                             "poplar_SNPs_Chr10_14460to14550kb.traw",
                             package = "SKATMCMT")

# for(raw_file_path in list.files(
#   "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/00_SNP_format_conversions/wholeChr1323geno/",
#   full.names = TRUE)){
#   print(paste("Run MTMC-SKAT for scaffold", raw_file_path))
#   MTMCSKAT_workflow(phenodata = phenodata,
#                     #phenodata = "/scratch2/NSF_GWAS/SKAT_SLURMS/inputs/shoot_5w.header.pheno",
#                     covariates = covariates,
#                     raw_file_path = raw_file_path,
#                     window_size = 3000,
#                     window_shift = 1000,
#                     output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
#                     pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/mtmcskat/pre_allocated_dir",
#                     n_core = "AllCores",
#                     #n_core = 2,
#                     job_id = "my_sample_analysis_110320_remove_callSKAT",
#                     desired_sig_figs = 2,
#                     min_accuracy = 2,
#                     max_accuracy = 5,
#                     plot = TRUE)
# }

MTMCSKAT_workflow(phenodata = phenodata,
                  #phenodata = "/scratch2/NSF_GWAS/SKAT_SLURMS/inputs/shoot_5w.header.pheno",
                  covariates = covariates,
                  raw_file_path = raw_file_path,
                  window_size = 3000,
                  window_shift = 1000,
                  output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
                  pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/mtmcskat/pre_allocated_dir",
                  n_core = "AllCores",
                  #n_core = 2,
                  job_id = "my_sample_analysis_110320_remove_callSKAT",
                  desired_sig_figs = 2,
                  min_accuracy = 2,
                  max_accuracy = 5,
                  plot = TRUE)
                  #RAM = 4e9)

