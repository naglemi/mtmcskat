#' Multithreaded calculation of null p-values to enable calculation of empirical
#' p-values
#'
#' These functions are used to calculate empirical p-values given initial
#' p-values and null p-values from `SKAT`.
#'
#' Prior to completion of jobs by cores in `MTMC-SKAT`, the number of p-values
#' above or below the in initial p-value are calculated using `tally_p_null`.
#' This reduces the load of data to be returned. Empirical p-values are then
#' calculated from all returned data using `p_empirical_from_tally`.
#'
#' @param p_table a data frame or matrix with a row for each SNP window and a
#'   column for each p-value resulting from a permuted null model
#'
#' @return a data frame with four columns: 1) position, in base pairs, for
#' the center of a SNP; 2) the original p-value
#'   for the given SNP window, obtained by SKAT without resampling; 3) The
#'   number of permuted p-values from a given thread that are above the
#'   original p-value for a given window; and 4) The number of permuted
#'   p-values from a given thread that are below the original p-value for a
#'   given window
#' @export
#'
#' @examples
#' \dontrun{sample_p_null_tallies <- tally_p_null(sample_p_table)}
#'
tally_p_null <- function(p_table){

  p_table <- stats::na.omit(p_table)

  empirical_p_table <- as.data.frame(cbind(p_table[,1:2],
                                           rep(NA, nrow(p_table)),
                                           rep(NA, nrow(p_table))))

  colnames(empirical_p_table) <- c("position", "SKAT_p-val",
                                   "n_perm_above", "n_perm_ltoreq")

  for(k in 1:nrow(p_table)){

    n_perm_ltoreq <- length(which(
      p_table[k, 3:ncol(p_table)] <= unlist(p_table[k, 2])))

    n_perm_above <- length(which(
      p_table[k, 3:ncol(p_table)] > unlist(p_table[k , 2])))

    empirical_p_table$n_perm_ltoreq[k] <- n_perm_ltoreq
    empirical_p_table$n_perm_above[k] <- n_perm_above
  }

  empirical_p_table
}

#' Tally p-values from permuted null model to enable calculation of empirical
#' p-values
#'
#' @param p_null_tallies a matrix containing tallies of null p-values above or
#'   below the initial p-value
#' @param scaffold_ID a string or numeric value indicating the scaffold, to be
#'   used only for labeling output
#' @param min_n_permutations Integer, the minimum number of permutations needed
#'   to calculate the empirical p-value out to the desired number of significant
#'   figures
#'
#' @return A dataframe with four columns, for 1) scaffold ID, 2) SNP window
#' position, 3) p-values from the model used in SKAT without resampling, and
#' 4) empirical p-values
#' @export
#'
#' @examples
#' data("sample_p_null_tallies")
#' p_empirical_from_tally(p_null_tallies = sample_p_null_tallies,
#'                        scaffold_ID = 10,
#'                        min_n_permutations = 1e5)
p_empirical_from_tally <- function(p_null_tallies, scaffold_ID,
                                   min_n_permutations){
  total_perm_p_ltoreq <- stats::aggregate(
    p_null_tallies$n_perm_ltoreq,
    by=list(p_null_tallies$position,
            p_null_tallies$`SKAT_p-val`),
    FUN=sum) # https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group

  total_perm_p_ltoreq$x <- total_perm_p_ltoreq$x

  total_perm_p_above <- stats::aggregate(
    p_null_tallies$n_perm_above,
    by=list(p_null_tallies$position,
            p_null_tallies$`SKAT_p-val`), FUN=sum)

  total_perm_p_ltoreq$empirical_p <- (total_perm_p_ltoreq$x + 1) /
    ( total_perm_p_ltoreq$x + total_perm_p_above$x + 1)

  output <- cbind(scaffold_ID,
                  total_perm_p_ltoreq$Group.1,
                  total_perm_p_ltoreq$Group.2,
                  total_perm_p_ltoreq$empirical_p)

  output <- as.data.frame(output)

  colnames(output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled")

  output$`SKAT_p-val` <- as.numeric(as.character(output$`SKAT_p-val`))
  output$`SKAT_p-val_resampled` <- as.numeric(as.character(
    output$`SKAT_p-val_resampled`))

  # round to significant figures (which are limtied by n_permutations)
  output$`SKAT_p-val_resampled` <- round(
    output$`SKAT_p-val_resampled`,
    digits = log(min_n_permutations,
                 base = 10))

  output
}

