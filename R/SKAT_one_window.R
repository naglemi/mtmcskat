#' Title
#'
#' @param pos_and_SNPs
#' @param this_position
#' @param window_size
#' @param Z
#' @param raw_file_path
#' @param null_model
#' @param n_permutations
#' @param SKAT_O
#' @param resampling
#' @param return_all_p_vals
#'
#' @return
#' @export
#'
#' @importFrom SKAT SKAT
#'
#' @examples
SKAT_one_window <- function(this_position,
                            window_size,
                            Z,
                            raw_file_path,
                            null_model,
                            n_permutations,
                            SKAT_O = "OFF",
                            resampling=FALSE,
                            return_all_p_vals=FALSE){

  # Only run SKAT over a SNP window if the window actually contains SNPs
  if(is.matrix(Z)==TRUE){

    # (Before calculating empirical p-val, if applicable,)
    # ...calculate p-value for this SNP window's effect given model assumptions
    this_SKAT_out <- SKAT::SKAT(Z, null_model)

    # resampling is performed for mtmcskat but not mtskat
    if(resampling==TRUE){

      # For mtmcskat parallelized over SNPs, we calculate the
      # empirical p-value here and return it, since this can be done on one
      # thread for each empirical p-value
      if(return_all_p_vals == FALSE){

        p_empirical <- calculate_SKAT_empirical_p(
          Z = Z,
          n_permutations = n_permutations,
          null_model = null_model)

        to_append <- c(raw_file_path,
                       this_position,
                       as.numeric(as.character(this_SKAT_out$p.value)),
                       p_empirical, NA, NA)
      }

      # for mtmcskat parallelized over null models, we do not have a single
      # thread containing all permuted p-values, thus we cannot calculate
      # the empirical p-value without first returning all p-values into
      # one place
      if(return_all_p_vals == TRUE){

        p_list <- calculate_SKAT_empirical_p(Z = Z,
                                             n_permutations = n_permutations,
                                             null_model = null_model,
                                             return_all_p = TRUE)

        to_append <- as.numeric(
          c(this_position,
            as.numeric(
              as.character(this_SKAT_out$p.value)), p_list))
      }

    }

    # no resampling for mtskat... we only need to generate initial p-values
    if(resampling==FALSE){

      p_empirical <- NA
      to_append <- c(raw_file_path,
                     this_position,
                     as.numeric(as.character(this_SKAT_out$p.value)),
                     p_empirical, NA, NA)
    }

  } else {
    message(paste0("This chunk's window at ",
                   this_position, " is not matrix!"))

    to_append <- as.data.frame(t(
      c(raw_file_path, this_position, NA, NA, NA, NA)))

  }

  to_append
}
