#' Select windows corresponding to p-values within range of interest
#'
#' @param data An output from \code{\link{mtskat}) or a similarly formatted
#'   data.frame
#' @param upper_bound A numeric value corresponding to the maximum p-value for
#'   the range of interest
#' @param lower_bound A numeric value correspondign to the minimum p-value for
#'   the range of interest
#'
#' @return A numeric vector with indices for the position (center of window) for
#'   each SNP window with p-values within the range of interest
#' @export
#'
#' @examples
#' select_windows_range_p(
#'   data = sample_mtskat_results,
#'   upper_bound = 0.01,
#'   lower_bound = 0.001)
#'
select_windows_range_p <- function(data,
                                   upper_bound,
                                   lower_bound){

  window_list <- data[which(
    data$`SKAT_p-val` < upper_bound & data$`SKAT_p-val` >= lower_bound),
    ]$position

  window_list
}
