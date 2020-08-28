#' Call SKAT
#' Merge this function in with SKAT_one_window
#'
#' Run SKAT with or without resampling, according to user specfiications. In the latter case, also calculate empirical p-value.
#'
#' @param pos_and_SNPs
#' @param scaffold_ID
#' @param resampling
#' @param null_model
#' @param n_permutations
#' @param return_all_p_vals
#'
#' @return
#' @export
#'
#' @examples
#' callSKAT(pos_and_SNPs = pos_and_SNPs[[i]], scaffold_ID = scaffold_ID, resampling = TRUE, null_model = null_model, n_permutations = n_permutations)
#'
callSKAT <- function(pos_and_SNPs,
                     scaffold_ID,
                     resampling,
                     null_model,
                     n_permutations,
                     return_all_p_vals = FALSE){

  # this_position <- pos_and_SNPs[[1]]
  # Z <- pos_and_SNPs[[2]]
  #browser()
  if(is.matrix(pos_and_SNPs[[2]])==TRUE){
    #print("About to go into SKAT_one_window")
    to_append <- SKAT_one_window(this_position = pos_and_SNPs[[1]],
                                 #window_size,
                                 Z = pos_and_SNPs[[2]],
                                 raw_file_path = scaffold_ID,
                                 resampling = resampling,
                                 null_model = null_model,
                                 n_permutations = n_permutations)
    #print("Out of SKAT_one_window")
  } else {
    message(paste0("This chunk's window at ", pos_and_SNPs[[1]], " is not matrix!"))
    to_append <- as.data.frame(t(c(scaffold_ID, pos_and_SNPs[[1]], NA, NA, NA, NA)))

    #browser()
  }

  #print(paste0("From within callSKAT, will append..."))
  #print(to_append)
  to_append
}

