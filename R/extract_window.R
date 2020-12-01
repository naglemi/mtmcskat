#' Extract and format a window of SNPs to be tested in a SKAT kernel
#'
#' @param this_position An integer, indicating the center of a SNP window for
#'   which the user wishes to extract SNPs (in base pairs)
#' @param window_size An integer, indicating the size of each SNP window (in
#'   base pairs)
#' @param genodata A matrix obtained by reading data from `.traw` format (see
#'   (see \url{https://www.cog-genomics.org/plink2/formats}{PLINK
#'   documentation})) into R
#'
#' @return A list containing three objects: an integer of chromosome of origin,
#'   an integer of SNP window center position, and a matrix of alternative
#'   allele counts for each genotype and SNP within the given window
#' @export
#'
#' @examples
#' small_genodata_path <- system.file("extdata",
#'                                    "poplar_200genotypes_14490to14520kb.traw",
#'                                    package = "mtmcskat")
#' small_genodata <- data.table::fread(small_genodata_path)
#'
#' extract_window(this_position = 145e5,
#'   window_size = 3000,
#'   genodata = small_genodata)
#'
extract_window <- function(this_position, window_size, genodata){

  locus_of_interest <- this_position
  window_start <- as.numeric(as.character(locus_of_interest)) - (window_size/2)

  if (window_start < 0){
    window_start <- 0
  }
  window_end <- as.numeric(as.character(locus_of_interest)) + (window_size/2)

  indices_to_pull <-
    which(
      with(genodata,
           genodata$POS >= window_start & genodata$POS <= window_end),
      arr.ind = TRUE)

  if(length(indices_to_pull) >= 1){
    genodata_thiswindow <-
      genodata[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)],]

    Chr <- unique(genodata_thiswindow[,1])
    if(length(Chr) < 1) {
      stop("Where is chromosome data?")
    }
    if(length(Chr) > 1) {
      stop(paste(Sys.time(),
                 " - extract_window should be provided with no more than a",
                 "single scaffold, but appears to have multiple."))
    }

    Z <- convert_to_Z(genodata_thiswindow = genodata_thiswindow)

    out_list <- list(this_position, Z, Chr)
  }

  if(length(indices_to_pull) == 0) {
    message(paste0("No SNPs within ",window_size/1000,"kb window with center",
                   " of ", this_position))
    out_list <- list(NA, NA, NA)
  }


  names(out_list) <- c("Position", "Z", "Chr")

  return(out_list)
}
