iter_window <- function(this_scaff_subset, times, window_list, window_size) {

  nextEl <- function() {
    if (times < length(window_list)) {
      times <<- times + 1
    }
    else {
      stop('StopIteration')
    }

    window_boundaries <- bound_window(this_position = window_list[times],
                                      window_size = window_size)

    # print(paste0("Window bounded from ",
    #              window_boundaries$window_start,
    #              " to ",
    #              window_boundaries$window_end))
    #print(window_list[times])
    genodata_thiswindow <- pull_window(this_scaff_subset = this_scaff_subset,
                                       window_start = window_boundaries$window_start,
                                       window_end = window_boundaries$window_end)

    if(is.null(genodata_thiswindow)) {
      Z <- NA
    } else {
      Z <- convert_to_Z(genodata_thiswindow = genodata_thiswindow)
    }

    pos_and_Z <- list(position=window_list[times], Z=Z)

    pos_and_Z
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('iter_window', 'abstractiter', 'iter')

  return(obj)
}
