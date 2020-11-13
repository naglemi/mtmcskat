#' Divide a list of SNP windows into a chunk for each thread
#'
#' @inheritParams re_allocate_windows
#' @param n_thread An integer indicating the number of threads to be used for
#'   multithreading
#'
#' @return a list of lists of lists, in which the list of lists passed by the
#'   `pre_allocated_SNP_windows` argument is divided into a number of lists
#'   indicated by the `n_threads` argument
#' @export
#'
#' @examples
#' data("sample_pre_allocated_SNP_windows")
#' pre_allocated_SNP_windows <- chunk_windows(
#'   pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows,
#'   n_thread = 2)
#'
chunk_windows <- function(pre_allocated_SNP_windows,
                          n_thread){

  chunk_size <- length(pre_allocated_SNP_windows) / n_thread

  message(paste0(Sys.time(), " - Chunking list of length ",
                 length(pre_allocated_SNP_windows),
                 " into blocks each with no more than",
                 ceiling(chunk_size),
                 " SNP windows each"))

  pre_allocated_SNP_windows <- split(
    pre_allocated_SNP_windows,
    (seq_along(pre_allocated_SNP_windows) - 1) %/% chunk_size)

  message(paste0("...resulting in ",
                 length(pre_allocated_SNP_windows),
                 " blocks"))

  pre_allocated_SNP_windows
}
