#' Benchmark the amount of RAM needed per permutation in a null model
#'
#' @param phenotype A phenotype file read into R from standard PLINK phenotype
#'   format (see (\href{https://www.cog-genomics.org/plink2/formats}{PLINK
#'   documentation}), with labeled columns for FID, IID, and trait.
#' @param covariates A covariate file read into R from either standard PLINK
#'   covariate format ((\href{https://www.cog-genomics.org/plink2/formats}{PLINK
#'   documentation}) or in a `.csv` with headers for genotype ID (in the first
#'   column) and an ID for each covariate (subsequent columns)
#' @param benchmark_size The number of permutations to be used for the benchmark
#'
#' @return Integer, the amount of RAM needed per permutation in a null model (in
#'   bytes)
#' @export
#'
#' @examples
#' data("sample_phenotype")
#' data("sample_covariates")
#'
#' nm_RAM_per_perm <- benchmark_nm(phenotype = sample_phenotype,
#'                                 covariates = sample_covariates,
#'                                 benchmark_size = 1000)
#'
benchmark_nm <- function(phenotype, covariates, benchmark_size = 1000){

  null_model <- SKAT::SKAT_Null_Model(phenotype ~ 1 + as.matrix(covariates),
                                      out_type="C",
                                      n.Resampling = benchmark_size,
                                      type.Resampling="bootstrap")

  nm_RAM_per_perm <- utils::object.size(null_model) / benchmark_size

  nm_RAM_per_perm
}

#' Calculate the number of permutations in a null model that can fit in RAM
#' partitioned for a single thread
#'
#' @param nm_RAM_per_perm Integer, the amount of RAM needed per permutation in a
#'   null model (in bytes); this is the output from \code{\link{benchmark_nm}}
#' @param RAM Integer, the total amount of RAM available to be used, for all
#'   threads (in bytes)
#' @param n_thread Integer, the number of threads over which RAM will be divided
#' @param verbose If TRUE, print messages explaining logic
#'
#' @return An integer indicating the amount of RAM that is available to be used
#'   by any given thread
#' @export
#'
#' @examples
#' data("sample_phenotype")
#' data("sample_covariates")
#'
#' nm_RAM_per_perm <- mtmcskat::benchmark_nm(phenotype = sample_phenotype,
#'                                           covariates = sample_covariates,
#'                                           benchmark_size = 1000)
#'
#' RAM <- benchmarkme::get_ram()
#'
#' calculate_max_perm_per_core(nm_RAM_per_perm = nm_RAM_per_perm,
#'                             RAM = RAM,
#'                             n_thread = 2)
#'
calculate_max_perm_per_core <- function(nm_RAM_per_perm,
                                        RAM,
                                        n_thread,
                                        verbose = TRUE){

  # What is the maximum number of permutations that can be performed without a
  #   round of boss-worker communication?
  max_simultaneous_perm <- floor(as.numeric(RAM / nm_RAM_per_perm))

  max_simultaneous_perm_per_core <- floor(max_simultaneous_perm / n_thread)

  if(verbose == TRUE){
    message(paste("Amount of RAM we are banking on is", RAM/1e9,
                  "GB and each permutation requires", nm_RAM_per_perm/1e3, "kB,",
                  "so there can be",
                  format(max_simultaneous_perm,
                         big.mark=",",scientific=FALSE),
                  "permutations at once. Divided over", n_thread,
                  "threads, each thread can run up to",
                  format(max_simultaneous_perm_per_core,
                         big.mark=",",scientific=FALSE), "permutations.\n"))
  }

  max_simultaneous_perm_per_core
}

size_RAM_wiggle <- function(RAM, wiggle_factor){

  RAM <- RAM/wiggle_factor
  RAM
}

