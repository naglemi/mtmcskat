#' Obtain SKAT results, including permuted or empirical p-values, for a SNP
#' window or chunk of SNP windows
#'
#' @param pos_and_SNPs A list containing SNP data divided into windows, along
#'   with metadata (scaffold and SNP window center position) produced by
#'   \code{\link{pre_allocate}}, and if passing chunks of windows rather than
#'   individual windows to each worker, broken into chunks using
#'   \code{\link{chunk_windows}}.
#' @param chunk Should be TRUE if the list passed by `pos_and_SNPs` argument has
#'   been processed into chunks by `chunk_windows()`, as is done in a standard
#'   \code{\link{mtskat}} workflow
#' @inheritParams SKAT_one_window
#'
#' @return A dataframe to be appended to results. If chunking (see above), this
#'   will contain multiple rows, one for each SNP window. Otherwise, it will
#'   contain a single row for one SNP window. Each row will contain a vector
#'   produced by \code{\link{SKAT_one_window}} (see documentation for that
#'   function for further details).
#' @export
#'
#' @examples
#' mappable_SKAT(
#'   pos_and_SNPs = sample_pre_allocated_SNP_windows[[1]],
#'   scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]],
#'   null_model = sample_null_model,
#'   resampling = TRUE,
#'   n_permutations = 1000,
#'   chunk = FALSE)
#'
#'
mappable_SKAT <- function(pos_and_SNPs,
                          scaffold_ID,
                          null_model,
                          n_permutations,
                          resampling=FALSE,
                          return_all_p_vals=FALSE,
                          chunk = TRUE){

  if ( chunk == TRUE ) {

    result_df <- as.data.frame(matrix(NA, nrow=1, ncol=4))

    for ( i in 1:length(pos_and_SNPs)){

      to_append <- SKAT_one_window(this_position = pos_and_SNPs[[i]][[1]],
                                   Z = pos_and_SNPs[[i]][[2]],
                                   scaffold_ID = scaffold_ID,
                                   resampling = resampling,
                                   null_model = null_model,
                                   n_permutations = n_permutations)

      result_df <- rbind(result_df, to_append)

    }
    return(result_df[-1,])
  }

  if ( chunk == FALSE ) {
    to_append <- SKAT_one_window(this_position = pos_and_SNPs[[1]],
                                 Z = pos_and_SNPs[[2]],
                                 scaffold_ID = scaffold_ID,
                                 resampling = resampling,
                                 null_model = null_model,
                                 n_permutations = n_permutations)

    #pos_and_SNPs <-NULL
    #gc()
    return(as.data.frame(t(to_append)))

  }
}
