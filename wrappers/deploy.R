# SKAT_mydata.R
#rm(list=ls())

#install.packages("SKAT")
#install.packages("tari")
library(SKAT)
library(data.table)
library(tools)
library(foreach)
library(doParallel)
library(plyr)

library(optparse)
library(iterators)
library(itertools)
library(furrr)
library(purrr)

library(dplyr)
library(qqman)
library(ggplot2)

library(SKATMCMT)
library(grid)
library(gridGraphics)

#library(SKATMCMT)

# Read arguments from command line  ---------------------------------------


option_list = list(
  make_option(c("-p", "--phenodata_path"),
              type="character",
              default=NULL,
              help="path to phenotype file",
              metavar="character"),
  make_option(c("-c", "--covariates"),
              type="character",
              default=NULL,
              help="path to covariates file",
              metavar="character"),
  make_option(c("-r", "--raw_file_path"),
              type="character",
              default=NULL,
              help="path to genodata in raw format (see PLINK manual)",
              metavar="character"),
  make_option(c("-S", "--window_size"),
              type="numeric",
              default=3000,
              #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
              metavar="numeric"),
  make_option(c("-s", "--window_shift"),
              type="numeric",
              default=1000,
              #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
              metavar="numeric"),
  make_option(c("-o", "--output_dir"),
              type="character",
              default="/scratch2/NSF_GWAS/Results/SKAT/",
              help="Path to directory where results will be saved",
              metavar="character"),
  make_option(c("-j", "--job_id"),
              type="character",
              default="no_job_id",
              help="User-provided job identifier, to be used as subfolder in which results are placed",
              metavar="character"),
  make_option(c("-n", "--ncore"),
              type="character",
              default="AllCores",
              help="A numeric value representing the number of cores to be used, or a string AllCores to automatically detect and use all cores",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# opt <- list(phenodata_path = "/scratch2/NSF_GWAS/phenodata/third_training/PCs_ratios_binarized_4class/TDZ_shoot_area.emmax.pheno",
#             #phenodata_path = "/scratch2/NSF_GWAS/phenodata/epicormic/EpicormicBudBreak_South_March_2014.header.pheno",
#             covariates = "/scratch2/NSF_GWAS/phenodata/Covariates_DiamFinal_PhaseFinal_PCs.txt",
#             #raw_file_path = "/scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr10.traw",
#             raw_file_path = "/scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr10_portion.traw",
#             #raw_file_path = "/scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr10_10kb.traw",
#             window_size = 3000,
#             window_shift = 1000,
#             output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
#             job_id = "fixing_iter_a42",
#             ncore = 48)

setwd(opt$output_dir)

phenodata <- fread(opt$phenodata_path)
this_phenotype <- as.numeric(unlist(phenodata[,3]))

covariates <- fread(opt$covariates)

whole_genome_start <- proc.time()


# cl <- makeCluster(availableCores(),
#                   type = "FORK",
#                   outfile="cluster.outfile")
# registerDoParallel(cl)
# on.exit(stopCluster(cl))
#plan(cl)
plan('multicore')

RAM_GB <- 5
options(future.globals.maxSize= RAM_GB*1000*(1024^2)) # https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r

runSKATtraw(phenodata = phenodata,
            covariates = covariates,
            raw_file_path = opt$raw_file_path,
            window_size = opt$window_size,
            window_shift = opt$window_shift,
            output_dir = opt$output_dir,
            job_id = opt$job_id,
            chunk_size = 500,
            ncore = opt$ncore)



# Stop the clock
proc.time() - whole_genome_start

