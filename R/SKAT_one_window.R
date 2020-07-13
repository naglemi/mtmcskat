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
SKAT_one_window <- function(pos_and_SNPs, this_position, window_size, Z, raw_file_path, null_model, n_permutations, SKAT_O = "OFF", resampling=FALSE, return_all_p_vals=FALSE){
  ##print("Inside SKAT_one_window")
  # if(hasArg(pos_and_SNPs)){
  #   this_position <- pos_and_SNPs[[1]]
  #   Z <- pos_and_SNPs[[2]]
  # }
  #browser()
  this_SKAT_out <- SKAT(Z, null_model, kernel = "linear.weighted")
  ##print("Done running SKAT")
  if(resampling==TRUE){
    if(return_all_p_vals == FALSE){
      p_empirical <- calculate_SKAT_empirical_p(Z = Z,
                                                n_permutations = n_permutations,
                                                null_model = null_model)
      to_append <- c(raw_file_path, this_position, as.numeric(as.character(this_SKAT_out$p.value)), p_empirical, NA, NA)
    }
    if(return_all_p_vals == TRUE){
      #print("Running SKAT to generate p-list")
      p_list <- calculate_SKAT_empirical_p(Z = Z,
                                           n_permutations = n_permutations,
                                           null_model = null_model,
                                           return_all_p = TRUE)
      #browser()
      #print("Generating output to append")
      to_append <- as.numeric(c(this_position, as.numeric(as.character(this_SKAT_out$p.value)), p_list))
    }

  }
  if(resampling==FALSE){
    ##print("Not resampling (yet)")
    p_empirical <- -9
    to_append <- c(raw_file_path, this_position, as.numeric(as.character(this_SKAT_out$p.value)), p_empirical, NA, NA)
  }

  #print(length(to_append))
  ##print(to_append)

  to_append
}
