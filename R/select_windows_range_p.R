#' Select windows corresponding to p-values within range of interest
#'
#' @param x An output from \code{\link{mtskat}} or a similarly formatted
#'   x.frame
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
#' data("small_mtskat_results")
#'
#' select_windows_range_p(
#'   x = small_mtskat_results,
#'   upper_bound = 0.01,
#'   lower_bound = 0.001)
#'
select_windows_range_p <- function(x,
                                   upper_bound,
                                   lower_bound){

  # Most (if not all) windows tested in each round are being tested because
  #   their initial model p-value, from SKAT without resampling, indicates that
  #   the specific # permutations we're now working with is the right number
  #   needed to calculate their empirical p-values with desired accuracy
  window_list <- x[which(
    x$`SKAT_p-val` < upper_bound & x$`SKAT_p-val` >= lower_bound),
    ]$position

  # For what windows does the empirical p-value have more leading zeros than
  #   expected, due to our "guess" of the # permutations needed, based on the
  #   initial model p-value, being too low? We will test these again with more
  #   permutations. Note, I have never observed this theoretical case
  #   actually happening, but it is in tests to make sure functionality works.
  window_list_2 <- x[which(
    x$`SKAT_p-val_resampled` < upper_bound &
      x$`SKAT_p-val_resampled` >= lower_bound),
    ]$position

  all_windows <- c(window_list, window_list_2)

  if(length(all_windows) != length(unique(all_windows))){
    stop("Error: Duplicate windows in `select_windows_range_p")
  }

  all_windows
}
