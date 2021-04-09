#' Pre-allocate a list of overlapping SNP windows
#'
#' SNP data is converted into overlapping windows as specified. This data
#' structure preparation is useful for parallelization of
#' \code{\link[SKAT]{SKAT}}.
#'
#' @param raw_file_path complete file path to SNP data in `.traw` format (see
#'   \href{https://www.cog-genomics.org/plink2/formats}{PLINK documentation})
#' @param window_shift An integer, indicating the number of base pairs over
#'   which each rolling window will slide; in other terms, the distance between
#'   the start (or end) positions of adjacent overlapping windows
#' @param pre_allocated_dir a directory where pre-allocated SNP window lists are
#'   kept
#' @inheritParams extract_window
#'
#' @return a list of lists, with each sub-list containing elements as described
#'   in documentation for \code{\link{extract_window}}
#' @export
#'
#' @examples
#'
#' \dontrun{
#' raw_file_path <- system.file("extdata",
#'   "poplar_SNPs_Chr10_14460to14550kb.traw",
#'   package = "SKATMCMT")
#'
#' pre_allocate(pre_allocated_dir = tempdir(),
#'   raw_file_path = raw_file_path,
#'   window_size = 3000,
#'   window_shift = 1000)
#' }
#'
pre_allocate <- function(raw_file_path, window_size, window_shift,
                         pre_allocated_dir,
                         impute_to_mean = TRUE,
                         remove_novar_SNPs = TRUE,
                         missing_cutoff = 0.15){

  if(!dir.exists(pre_allocated_dir)) dir.create(pre_allocated_dir,
                                                recursive = TRUE)

  pre_allocated_path <- paste0(
    pre_allocated_dir, "/",
    basename(tools::file_path_sans_ext(raw_file_path)), ".",
    window_size, "bp_win.rds")

  if(file.exists(pre_allocated_path)){
    message(paste0("Data structure already exists. Read in from: ",
                   pre_allocated_path,
                   "\n"))
    pos_and_SNP_list <- readRDS(pre_allocated_path)
  }
  if(!file.exists(pre_allocated_path)){
    genodata <- data.table::fread(paste0(raw_file_path),
                                           fill = TRUE,
                                  tmpdir = pre_allocated_dir)

    message(paste("Read genotype data for raw file", raw_file_path,
                  "with", dim(genodata)[1], "SNPs and",
                   dim(genodata)[2] - 6, "genotypes"))
    genodata$POS <- as.numeric(as.character(genodata$POS))

    max_window <- max(genodata$POS)

    # The starting position is either A) one half of the window size,
    #   since we declare each window position to be the center of each window
    #   and wish to start at the start of the genome; or
    #   B) the first SNP position rounded to one-half of the window size,
    #   for situations in which the first window declared by approach A
    #   would not contain any SNP.
    starting_position <- max((window_size / 2),
                             plyr::round_any(
                               min(genodata$POS), window_size / 2))

    window_list <- seq(starting_position,
                       max_window,
                       by=window_shift)

    message("Allocating data structure")
    ptm <- proc.time()

    pos_and_SNP_list <- lapply(window_list,
                               function(x) extract_window(
                                 this_position = x,
                                 window_size = window_size,
                                 genodata = genodata,
                                 impute_to_mean = impute_to_mean,
                                 remove_novar_SNPs = remove_novar_SNPs,
                                 missing_cutoff = missing_cutoff))

    time <- proc.time() - ptm
    message(paste("Took", time[3],
                  "seconds to allocate pos, SNP data structure"))

    saveRDS(pos_and_SNP_list, pre_allocated_path)

    message(paste("Pre-allocated SNP data saved to",
                  pre_allocated_path,
                  "\n"))
  }

  genodata <- NULL # Will this prevent it from being exported on slurm?
  gc()

  print(paste("Pre-allocated SNP window data takes up",
              utils::object.size(pos_and_SNP_list)/1e6,
              "MB\n"))

  pos_and_SNP_list
}
