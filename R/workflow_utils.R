#' Determine boundaries for subsetting SNP windows based on intial p-values
#'
#' @param leading_0s An integer for the number of 0's between the decimal point
#' and the first nonzero value in the maximum p-value for SNP windows that will
#' be included in the subset of interest
#' @param desired_sig_figs An integer for the number of significant figures
#' the user desires for empirical p-values
#' @param terminal_resampling If `TRUE`, the lower boundary will be 0 so all
#' remaining SNP windows will be included in the final set for resampling, even
#' if their empirical p-values cannot be calculated out to the user-specified
#' number of significant figures
#'
#' @return A list with two decimal values, named `upper` and `lower`
#' @export
#'
#' @examples
#'
#' # Produce boundaries from 0.01 to 0.001
#' determine_subset_pval_boundaries(leading_0s = 2,
#'                                  desired_sig_figs = 2,
#'                                  terminal_resampling = FALSE)
#'
#' # Produce boundaries from 0.01 to 0
#' determine_subset_pval_boundaries(leading_0s = 2,
#'                                  desired_sig_figs = 2,
#'                                  terminal_resampling = TRUE)
#'
#' # Produce boundaries from 0.01 to 0.0001
#' determine_subset_pval_boundaries(leading_0s = 2,
#'                                  desired_sig_figs = 3,
#'                                  terminal_resampling = FALSE)
#'
determine_subset_pval_boundaries <- function(leading_0s,
                                             desired_sig_figs,
                                             terminal_resampling){

  if(terminal_resampling == TRUE){
    lower_bound_p_val_for_MC_subset <- 0
  }
  if(terminal_resampling == FALSE){
    lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
  }
  upper_bound_p_val_for_MC_subset <- 0.1^leading_0s

  boundaries <- list("upper" = upper_bound_p_val_for_MC_subset,
                     "lower" = lower_bound_p_val_for_MC_subset)

  boundaries
}

#' Estimate, based on initial p-values, the number of permut
#'
#' This function estimates, based on the number of leading 0's in an initial
#' p-value, how many permutations will be needed to calculate an empirical
#' p-value out to the user-specified number of significant figures
#'
#' @inheritParams determine_subset_pval_boundaries
#'
#' @return An integer, representing the number of permutations expected to be
#' needed to calculate an empirical p-value with the desired level of accuracy
#' @export
#'
#' @examples
#' # 10,000 permutations should be needed to calculate p-values below 0.01
#' # (e.g. 0.0095 has two leading 0's) out to two significant figures
#' determine_n_permutations(leading_0s = 2, desired_sig_figs = 2)
determine_n_permutations <- function(leading_0s,
                                     desired_sig_figs){
  n_permutations <- 10^( leading_0s + desired_sig_figs )

  n_permutations
}
