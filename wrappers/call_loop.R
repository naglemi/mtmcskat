phenodata <- system.file("extdata",
                         "TDZ_shoot_area.plink.pheno",
                         package = "mtmcskat")

# phenodata <- "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/shoot_5w.header.pheno"
#
# covariates <- "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt"

covariates <- system.file("extdata",
                          "poplar_PCs_covariates.tbt",
                          package = "mtmcskat")

raw_file_path <- system.file("extdata",
                             "poplar_SNPs_Chr10_14460to14550kb.traw",
                             package = "mtmcskat")

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
                  covariates = covariates,
                  raw_file_path = raw_file_path,
                  window_size = 3000,
                  window_shift = 1000,
                  output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
                  pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/mtmcskat/pre_allocated_dir",
                  n_thread = "AllCores",
                  job_id = "my_sample_analysis_112920_mtmcskat_unit_tests",
                  desired_sig_figs = 2,
                  min_accuracy = 2,
                  max_accuracy = 5,
                  plot = TRUE)
                  #RAM = 4e9)

MTMCSKAT_workflow(phenodata = "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/shoot_5w.header.pheno",
                  covariates = "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt",
                  raw_file_path = "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/00_SNP_format_conversions/1323_cohort_mincount1_defaultmissingrates_Chr12.snp.pass.traw",
                  window_size = 3000,
                  window_shift = 1000,
                  output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
                  pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/pre_allocated_dir",
                  n_thread = 24,
                  job_id = "my_sample_analysis_111320_debug_moreefficient_NullModelpar",
                  desired_sig_figs = 2,
                  min_accuracy = 2,
                  max_accuracy = 5,
                  plot = TRUE,
                  RAM = 128e9)

MTMCSKAT_workflow(phenodata = "/oasis/projects/nsf/osu123/naglemi/test_one_phenotype/shoot_5w.header.pheno",
                  covariates = "/oasis/projects/nsf/osu123/naglemi/covariates/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt",
                  raw_file_path = "/oasis/projects/nsf/osu123/naglemi/test_one_scaffold_wholeChr_1323geno/1323_cohort_mincount1_defaultmissingrates_Chr12.snp.pass.traw.gz",
                  window_size = 3000,
                  window_shift = 1000,
                  output_dir = "/oasis/projects/nsf/osu123/naglemi/mtmcskat_out",
                  pre_allocated_dir = "/oasis/projects/nsf/osu123/naglemi/pre_allocated",
                  n_thread = "AllCores",
                  job_id = "my_sample_analysis_111320_debug_a11",
                  desired_sig_figs = 2,
                  min_accuracy = 2,
                  max_accuracy = 5,
                  plot = TRUE)
