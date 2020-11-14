#' Obtain a new list of SNP windows, based on p-values from an earlier list
#'
#' @inheritParams select_windows_range_p
#' @param pre_allocated_SNP_windows Output from \code{\link{pre_allocate}},
#'   which is a list of lists, with each sub-list containing elements as
#'   described in documentation for \code{\link{extract_window}}
#'
#' @return Output is the same format as for \code{\link{pre_allocate}}: a list
#'   of lists, with each sub-list containing elements as described in
#'   documentation for \code{\link{extract_window}}
#' @export
#'
#' @examples
#' data("sample_pre_allocated_SNP_windows")
#'
#' re_allocate_windows(
#'   x = sample_mtskat_results,
#'   upper_bound = 0.01,
#'   lower_bound = 0.001,
#'   pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows)
#'
re_allocate_windows <- function(x,
                                upper_bound,
                                lower_bound,
                                pre_allocated_SNP_windows){

  find_leading_0s <- function(x){
    # Thank you David Arenburg for this great one-liner
    # https://stackoverflow.com/questions/35553244/count-leading-zeros-between-the-decimal-point-and-first-nonzero-digit
    attr(regexpr("(?<=\\.)0+", x, perl = TRUE), "match.length")
  }

  window_list <- select_windows_range_p(
    x = x,
    upper_bound = upper_bound,
    lower_bound = lower_bound)

  if(length(window_list) == 0){

    message(paste0(Sys.time(),
                   " - No SNPs with p-vals within the range bounded from ",
                   lower_bound, " to ",
                   upper_bound))

    return("None")
  }

  new_pre_allocated_SNP_windows <- subset_list_of_lists(
    list_of_lists = pre_allocated_SNP_windows,
    desired_list = window_list,
    subindex = 1)

  # message(paste0("Number of SNP windows outside of subset_list_of_lists is ",
  #                length(new_pre_allocated_SNP_windows)))

  if(length(new_pre_allocated_SNP_windows) == 0){
    stop("Error: There is no SNP data for windows in the list provided.")
  }

  new_pre_allocated_SNP_windows

}
