iter_window_de <- function(this_scaff_subset, this_position, max_window, window_list, window_size, window_shift) {

  nextEl <- function() {
    if (this_position <= max_window) {
      this_position <<- this_position + window_shift
    }
    else {
      stop('StopIteration')
    }

    window_boundaries <- bound_window(this_position = this_position,
                                      window_size = window_size)

    genodata_thiswindow <- pull_window(this_scaff_subset = this_scaff_subset,
                                       window_start = window_boundaries$window_start,
                                       window_end = window_boundaries$window_end)

    Z <- convert_to_Z(genodata_thiswindow = genodata_thiswindow)

    pos_and_Z <- list(position=this_position, Z=Z)

    pos_and_Z
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('iter_window', 'abstractiter', 'iter')

  return(obj)
}
