# This is part of my latest attempt to parallelize at a high level over
#   scaffolds and phenotypes on SLURMS. We first need a callable script (this),
#   which will be called for each phenotype/scaffold combination in a separate
#   SLURMS batch job script. This approach is described here,
#   under "R Submission"
#   https://vsoch.github.io/lessons/sherlock-jobs/

library(optparse)
library(mtmcskat)

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
  make_option(c("-r", "--raw_file_path"),
              type="character",
              default=NA,
              help=paste("Path to folder of SNP data, organized into",
                         "chromosomes or other scaffolds, in PLINK .traw",
                         "format"),
              metavar="character"),
  make_option(c("--window_size"),
              type="numeric",
              default=3000,
              help=paste("Size of each SNP window (in base pairs)"),
              metavar="numeric"),
  make_option(c("--window_shift"),
              type="numeric",
              default=1000,
              help=paste("Shift between overlapping SNP windows",
                         "(in base pairs)"),
              metavar="numeric"),
  make_option(c("--output_dir"),
              type="character",
              default=NA,
              help=paste("Path to folder where outputs will be saved"),
              metavar="character"),
  make_option(c("--pre_allocated_dir"),
              type="character",
              default=NA,
              help=paste("Path to folder pre-allocated SNP data is saved",
                         "or will be saved if it doesn't already exist"),
              metavar="character"),
  make_option(c("--job_id"),
              type="character",
              default=NA,
              help=paste("Identifier for this group of jobs, to be submitted"),
              metavar="character"),
  make_option(c("--desired_sig_figs"),
              type="numeric",
              default=2,
              help=paste("numeric, number of significant figures out to which
                         the user desires empirical p-values be calculated"),
              metavar="numeric"),
  make_option(c("--min_accuracy"),
              type="numeric",
              default=2,
              help=paste("numeric, threshold x to begin resampling for",
                         "SNP windows with p-values below 10^-x; default",
                         "value is 2"),
              metavar="numeric"),
  make_option(c("--max_accuracy"),
              type="numeric",
              default=5,
              help=paste("numeric, limit x to end resampling for SNP windows",
                         "with p-values below 10^-x, usually due to",
                         "computational cost"),
              metavar="numeric"),
  make_option(c("--plot"),
              type="numeric",
              default=1,
              help=paste("If 1, make plots. If 0, make no plot."),
              metavar="numeric"),
  # make_option(c("-t", "--time"),
  #             type="character",
  #             default="1:00:00",
  #             help=paste("If running on SLURMS, a time string formatted as",
  #                        "hh:mm:ss indicating the maximum time to allow jobs",
  #                        "to run prior to termination"),
  #             metavar="character"),
  # make_option(c("-m", "--mode"),
  #             type="character",
  #             default="sequential",
  #             help=paste('character string of "sequential" or "slurms"'),
  #             metavar="character"),
  make_option(c("--RAM"),
              type="character",
              default=1,
              help=paste('Either "AllRAM" or the amount of RAM (in GB) which the
                         user wishes to set as a maximum for usable RAM in this
                         job'),
              metavar="character"),
  make_option(c("--n_thread"),
              type="numeric",
              default=24,
              help=paste('Number of threads for this job'),
              metavar="numeric"),
  make_option(c("--top_N"),
              type="numeric",
              default=Inf,
              help=paste('Integer representing the number of top associations',
              'on which the user wishes to perform resampling. For example,',
              '   if this value is set to 5, any SNPs that',
              '   do not produce p-values among the lowest 5 will not be',
              'included in outputs'),
              metavar="numeric"),
  make_option(c("--missing_cutoff"),
              type="numeric",
              default=0.15,
              help=paste('A numeric threshold representing the minimum desired',
                         'missing rate; missing rate is defined for each SNP',
                         'as the proportion of genotypes missing data for the',
                         'given SNP. Imputation to mean is performed , either',
                         'by `pre_allocate` or `SKAT` itself,',
                         'for all remaining missing values'),
              metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

MTMCSKAT_workflow(
  phenodata = as.character(opt$phenodata),
  covariates = as.character(opt$covariates),
  raw_file_path = as.character(opt$raw_file_path),
  window_size = opt$window_size,
  window_shift = opt$window_shift,
  output_dir = as.character(opt$output_dir),
  pre_allocated_dir = as.character(opt$pre_allocated_dir),
  job_id = opt$job_id,
  desired_sig_figs = opt$desired_sig_figs,
  min_accuracy = opt$min_accuracy,
  max_accuracy = opt$max_accuracy,
  plot = opt$plot,
  RAM = opt$RAM,
  n_thread = opt$n_thread,
  top_N = opt$top_N,
  missing_cutoff = opt$missing_cutoff)
