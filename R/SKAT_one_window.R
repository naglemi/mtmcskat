#' Run SKAT, optionally obtaining permuted p-values or empirical p-value, and
#' organize results into a vector to be added to dataframe of results
#'
#' @param this_position Integer indicating the center of a given SNP window (in
#'   base pairs)
#' @param scaffold_ID Integer indicating the chromosome or scaffold of interest
#' @param resampling If TRUE, will call \code{\link{calculate_SKAT_empirical_p}}
#'   to calculate an empirical p-value or obtain a vector of permuted p-values.
#'   If FALSE, will pass NA in place of an empirical p-value
#' @inheritParams calculate_SKAT_empirical_p
#'
#' @return If not returning all permuted p-values (as indicated by
#'   `return_all_p_vals`, the output will be a vector of length 4, containing
#'   `scaffold_ID`, `this_position`, the model p-value, and finally the
#'   empirical p-value (if applicable; otherwise NA). If returning all p-values,
#'   the length of the output vector will be two (for SNP scaffold and window
#'   center position) plus the number of permutations
#' @export
#'
#' @examples
#' data("small_pre_allocated_windows")
#' sample_null_model <- SKAT::SKAT_Null_Model(
#'   small_phenodata ~ 1 + as.matrix(small_covariates), out_type="C",
#'   n.Resampling = 1000)
#'
#' SKAT_one_window(
#'   this_position = small_pre_allocated_windows[[1]][[1]],
#'   Z = small_pre_allocated_windows[[1]][[2]],
#'   scaffold_ID = small_pre_allocated_windows[[1]][[3]],
#'   n_permutations = 1000,
#'   null_model = sample_null_model,
#'   resampling = TRUE,
#'   return_all_p = FALSE)


SKAT_one_window <- function(this_position,
                            Z,
                            scaffold_ID,
                            null_model,
                            n_permutations,
                            resampling=FALSE,
                            return_all_p_vals=FALSE,
                            missing_cutoff=0.15,
                            output_dir = NA,
                            output_basename = NA,
                            ...){

  # if(resampling == TRUE){
  #   browser()
  # }

  # Only run SKAT over a SNP window if the window actually contains SNPs
  if(is.matrix(Z)==TRUE){

    # (Before calculating empirical p-val, if applicable,)
    # ...calculate p-value for this SNP window's effect given model assumptions
    this_SKAT_out <- SKAT::SKAT(Z, null_model, missing_cutoff = missing_cutoff)

    # resampling is performed for mtmcskat but not mtskat
    if(resampling==TRUE){

      # For mtmcskat parallelized over SNPs, we calculate the
      # empirical p-value here and return it, since this can be done on one
      # thread for each empirical p-value
      if(return_all_p_vals == FALSE){

        p_empirical <- calculate_SKAT_empirical_p(
          Z = Z,
          n_permutations = n_permutations,
          null_model = null_model,
          missing_cutoff = missing_cutoff,
          output_dir = output_dir,
          output_basename = output_basename,
          this_position = this_position,
          scaffold_ID = scaffold_ID,
          ...)

        to_append <- c(scaffold_ID,
                       this_position,
                       as.numeric(as.character(this_SKAT_out$p.value)),
                       p_empirical)
      }

      # for mtmcskat parallelized over null models, we do not have a single
      # thread containing all permuted p-values, thus we cannot calculate
      # the empirical p-value without first returning all p-values into
      # one place
      if(return_all_p_vals == TRUE){

        p_list <- calculate_SKAT_empirical_p(Z = Z,
                                             n_permutations = n_permutations,
                                             null_model = null_model,
                                             return_all_p_vals = TRUE,
                                             missing_cutoff = missing_cutoff,
                                             output_dir = output_dir,
                                             output_basename = output_basename,
                                             this_position = this_position,
                                             scaffold_ID = scaffold_ID,
                                             ...)

        to_append <- as.numeric(
          c(this_position,
            as.numeric(
              as.character(this_SKAT_out$p.value)), p_list))
      }

    }

    # no resampling for mtskat... we only need to generate initial p-values
    if(resampling==FALSE){

      p_empirical <- NA
      to_append <- c(scaffold_ID,
                     this_position,
                     as.numeric(as.character(this_SKAT_out$p.value)),
                     p_empirical)
    }

  } else {
    #message(paste0("This chunk's window at ",
    #this_position, " is not matrix!"))

    to_append <- as.data.frame(t(
      c(scaffold_ID, this_position, NA, NA)))

  }

  #browser()

  to_append
}
