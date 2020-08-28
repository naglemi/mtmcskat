#' Title
#' Subset a list of lists down to a list of lists containing only lists with the value at a certain position (with this values indicated in another list)
#'
#' @param list_of_lists
#' @param desired_list list of desired elements
#' @param subindex index of lists within list of lists, which will be checked to see if the list contains the desired element and should be subsetted
#'
#' @return
#' @export
#'
#' @examples pre_allocated_SNP_windows <- subset_list_of_lists(list_of_lists = pre_allocated_SNP_windows, desired_list = window_list, subindex = 1)
subset_list_of_lists <- function(list_of_lists,
                                 desired_list,
                                 subindex = 1){

  # https://stackoverflow.com/questions/18556703/select-a-subset-of-lists-from-list-of-lists

  subsetted_list_of_lists <- list_of_lists[which(sapply(list_of_lists, `[[`, subindex) %in% desired_list)]

  subsetted_list_of_lists
}
