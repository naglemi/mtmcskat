#' Extract and format a window of SNPs to be tested in a SKAT kernel
#'
#' @param this_position An integer, indicating the center of a SNP window for
#'   which the user wishes to extract SNPs (in base pairs)
#' @param window_size An integer, indicating the size of each SNP window (in
#'   base pairs)
#' @param genodata A matrix obtained by reading data from `.traw` format (see
#'   (see \url{https://www.cog-genomics.org/plink2/formats}{PLINK
#'   documentation})) into R
#' @param impute_to_mean If `TRUE`, NA values for each SNP are replaced with
#'   the mean alternative allele count for the given SNP
#' @param remove_novar_SNPs If `TRUE`, SNPs with no variation will be removed
#' @param missing_cutoff A numeric threshold representing the minimum desired
#'   missing rate; missing rate is defined for each SNP as the proportion
#'   of genotypes missing data for the given SNP. Imputation to mean is
#'   performed , either by `pre_allocate` or `SKAT` itself,
#'   for all remaining missing values
#'
#' @return A list containing three objects: an integer of chromosome of origin,
#'   an integer of SNP window center position, and a matrix of alternative
#'   allele counts for each genotype and SNP within the given window
#' @export
#'
#' @examples
#' small_genodata_path <- system.file("extdata",
#'                                    "poplar_200genotypes_14490to14520kb.traw",
#'                                    package = "mtmcskat")
#' small_genodata <- data.table::fread(small_genodata_path)
#'
#' extract_window(this_position = 145e5,
#'   window_size = 3000,
#'   genodata = small_genodata)
#'
extract_window <- function(this_position, window_size, genodata,
                           impute_to_mean = TRUE,
                           remove_novar_SNPs = TRUE,
                           missing_cutoff = 0.15){

  locus_of_interest <- this_position
  window_start <- as.numeric(as.character(locus_of_interest)) - (window_size/2)

  if (window_start < 0){
    window_start <- 0
  }
  window_end <- as.numeric(as.character(locus_of_interest)) + (window_size/2)

  indices_to_pull <-
    which(
      with(genodata,
           genodata$POS >= window_start & genodata$POS <= window_end),
      arr.ind = TRUE)

  if(length(indices_to_pull) >= 1){
    genodata_thiswindow <-
      genodata[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)], ]

    Chr <- unique(genodata_thiswindow[, 1])
    if(length(Chr) < 1) {
      stop("Where is chromosome data?")
    }
    if(length(Chr) > 1) {
      stop(paste(Sys.time(),
                 " - extract_window should be provided with no more than a",
                 "single scaffold, but appears to have multiple."))
    }

    Z <- convert_to_Z(genodata_thiswindow = genodata_thiswindow)

    if(missing_cutoff > 0){
      missing_count_per_SNP <- colSums(is.na(Z))
      n_genotypes <- nrow(Z)

      missing_proportion_per_SNP <-
        missing_count_per_SNP /
        n_genotypes

      SNPs_with_missing_rate_gt_threshold <-
        which(missing_proportion_per_SNP > missing_cutoff, arr.ind = TRUE)

      if(length(SNPs_with_missing_rate_gt_threshold > 0)){
        Z <- Z[, -SNPs_with_missing_rate_gt_threshold]
      }
    }


    if(is.vector(Z)){ # Patch for cases when only 1 SNP in window
      Z <- as.matrix(Z)
    }
    if(ncol(Z) > 0) {

      if(remove_novar_SNPs == TRUE){
        # Count number of factors for each SNP
        factors_per_SNP <- c()
        for(i in 1:ncol(Z)){
          factors_per_SNP <- c(factors_per_SNP,
                               length(levels(factor(Z[, i]))))
        }

        n_SNPs <- ncol(Z)
        factors_per_SNP_gt1 <- which(factors_per_SNP > 1)
        if(n_SNPs != length(factors_per_SNP_gt1)){
          Z <- Z[, factors_per_SNP_gt1]
        }

      }

      if(is.vector(Z)){ # Patch for cases when only 1 SNP in window
        Z <- as.matrix(Z)
      }
      if(ncol(Z) > 0){ # if STILL > 0 SNPs after removing those w no variance ^
        if(impute_to_mean==TRUE){
          for(i in 1:ncol(Z)){
            Z[is.na(Z[, i]), i] <- mean(Z[, i], na.rm = TRUE)
          }
        }

        out_list <- list(this_position, Z, Chr)
      } else {
        out_list <- list(NA, NA, NA)
      }

    } else{
      out_list <- list(NA, NA, NA)
    }
  }

  if(length(indices_to_pull) == 0) {
    cat(paste0("No SNPs within ",window_size/1000,"kb window with center",
               " of ", this_position, "\n"))
    out_list <- list(NA, NA, NA)
  }


  names(out_list) <- c("Position", "Z", "Chr")

  return(out_list)
}
