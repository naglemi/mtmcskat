library(foreach)

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
  make_option(c("-a", "--pre_allocated_dir"),
              type="character",
              default=NA,
              help=paste("Path to folder pre-allocated SNP data is saved",
                         "or will be saved if it doesn't already exist"),
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
              metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

raw_file_path_list <- list.files(opt$genodata)

phenotype_files <- list.files(opt$phenodata)

foreach(raw_file_path = raw_file_path_list,
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
                             min_accuracy = opt$min_accuracy,
                             max_accuracy = opt$max_accuracy,
                             plot = TRUE,
                             RAM = "AllRAM",
                             n_core = "AllCores")

          pars
          }

print("Parameters to run: ")
print(pars)

if(opt$mode == "slurm"){
  library(rslurm)
  sjob <- slurm_apply(MTMCSKAT_workflow,
                      pars,
                      jobname = opt$job_id,
                      nodes = nrow(pars),
                      cpus_per_node = 1,
                      submit = FALSE,
                      slurm_options = list(time = opt$time))
}
if(opt$mode == "sequential"){
  apply(X = pars,
        MARGIN = 1,
        FUN = MTMCSKAT_workflow)
}


