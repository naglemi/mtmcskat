library(mtmcskat)
library(foreach)
library(optparse)

option_list = list(
  make_option(c("-p", "--phenodata"),
              type="character",
              default=NA,
              help="Path to folder of phenotype data files in PLINK format",
              metavar="character"),
  make_option(c("-c", "--covariates"),
              type="character",
              default=NA,
              help=paste("Path to a covariate file, comma-delimited with a",
                         "header, with genotypes in same order as phenotype",
                         "file"),
              metavar="character"),
  make_option(c("-g", "--genodata"),
              type="character",
              default=NA,
              help=paste("Path to folder of SNP data, organized into",
                         "chromosomes or other scaffolds, in PLINK .traw",
                         "format"),
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default=NA,
              help=paste("Path to folder where outputs will be saved"),
              metavar="character"),
  make_option(c("-l", "--pre_allocated_dir"),
              type="character",
              default=NA,
              help=paste("Path to folder pre-allocated SNP data is saved",
                         "or will be saved if it doesn't already exist"),
              metavar="character"),
  make_option(c("-j", "--job_id"),
              type="character",
              default=NA,
              help=paste("Identifier for this group of jobs, to be submitted"),
              metavar="character"),
  make_option(c("-S", "--window_size"),
              type="numeric",
              default=3000,
              help=paste("Size of each SNP window (in base pairs)"),
              metavar="numeric"),
  make_option(c("-s", "--window_shift"),
              type="numeric",
              default=1000,
              help=paste("Shift between overlapping SNP windows",
                         "(in base pairs)"),
              metavar="numeric"),
  make_option(c("-a", "--min_accuracy"),
              type="numeric",
              default=2,
              help=paste("numeric, threshold x to begin resampling for",
                         "SNP windows with p-values below 10^-x; default",
                         "value is 2"),
              metavar="numeric"),
  make_option(c("-A", "--max_accuracy"),
              type="numeric",
              default=5,
              help=paste("numeric, limit x to end resampling for SNP windows",
                         "with p-values below 10^-x, usually due to",
                         "computational cost"),
              metavar="numeric"),
  make_option(c("-t", "--time"),
              type="character",
              default="1:00:00",
              help=paste("If running on SLURMS, a time string formatted as",
                         "hh:mm:ss indicating the maximum time to allow jobs",
                         "to run prior to termination"),
              metavar="character"),
  make_option(c("-m", "--mode"),
              type="character",
              default="sequential",
              help=paste('character string of "sequential" or "slurms"'),
              metavar="character"),
  make_option(c("-n", "--n_thread"),
              type="numeric",
              default=24,
              help=paste('Number of threads for this job'),
              metavar="numeric")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# setwd("/scratch2/NSF_GWAS/notebooks/MTMC-SKAT")
# opt = list(phenodata = "test_two_phenotypes",
#            covariates = "/scratch2/NSF_GWAS/notebooks/InPlantaGWAS/02_Parsing_phenodata/pheno_files/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt",
#            genodata = "test_two_scaffolds_wholeChr_1323geno",
#            output = "/scratch2/NSF_GWAS/Results/SKAT/",
#            pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/mtmcskat/pre_allocated_dir",
#            job_id = "test_twogeno_twopheno_111020",
#            window_size = 3000,
#            window_shift = 1000,
#            min_accuracy = 2,
#            max_accuracy = 5,
#            mode = "sequential")

# opt = list(phenodata = "/oasis/projects/nsf/osu123/naglemi/test_two_phenotypes",
#            covariates = "/oasis/projects/nsf/osu123/naglemi/covariates/Covariates_AllPhasesButLast_Stem_noheader_PLINKformat.txt",
#            genodata = "/oasis/projects/nsf/osu123/naglemi/test_two_scaffolds_wholeChr_1323geno",
#            output = "/oasis/projects/nsf/osu123/naglemi/mtmcskat_out",
#            pre_allocated_dir = "/oasis/projects/nsf/osu123/naglemi/pre_allocated",
#            job_id = "firstrun_twoscaffolds_two_phenos_on_debug_node_a1",
#            window_size = 3000,
#            window_shift = 1000,
#            min_accuracy = 2,
#            max_accuracy = 5,
#            mode = "sequential")

raw_file_path_list <- list.files(opt$genodata,
                                 full.names = TRUE)

phenotype_files <- list.files(opt$phenodata,
                              full.names = TRUE)

pars <- foreach(raw_file_path = raw_file_path_list,
                .combine = rbind) %:%
  foreach(phenotype = phenotype_files,
          .combine = rbind) %do% {

          pars <- data.frame(phenodata = phenotype,
                             covariates = opt$covariates,
                             raw_file_path = raw_file_path,
                             window_size = opt$window_size,
                             window_shift = opt$window_shift,
                             output_dir = opt$output,
                             pre_allocated_dir = opt$pre_allocated_dir,
                             job_id = opt$job_id,
                             desired_sig_figs = 2,
                             min_accuracy = opt$min_accuracy,
                             max_accuracy = opt$max_accuracy,
                             plot = TRUE,
                             RAM = "AllRAM",
                             n_thread = opt$n_thread)
          }

if(opt$mode == "slurm"){
  library(rslurm)
  sjob <- slurm_apply(MTMCSKAT_workflow,
                      pars,
                      jobname = opt$job_id,
                      nodes = nrow(pars),
                      cpus_per_node = 1,
                      submit = FALSE,
                      slurm_options = list(time = opt$time,
                                           ntasks-per-node = opt$n_thread))
}

if(opt$mode == "sequential"){
  mapply(FUN = MTMCSKAT_workflow,
         phenodata = as.character(pars$phenodata),
         covariates = as.character(pars$covariates),
         raw_file_path = as.character(pars$raw_file_path),
         window_size = pars$window_size,
         window_shift = pars$window_shift,
         output_dir = as.character(pars$output_dir),
         pre_allocated_dir = as.character(pars$pre_allocated_dir),
         job_id = pars$job_id,
         desired_sig_figs = pars$desired_sig_figs,
         min_accuracy = pars$min_accuracy,
         max_accuracy = pars$max_accuracy,
         plot = pars$plot,
         RAM = pars$RAM,
         n_thread = "AllCores")
}


