#' Pre-allocate a list of overlapping SNP windows
#' SNP data is converted into overlapping windows as specified. This data structure preparation is useful for parallelization of (\code{"SKAT"}).
#'
#' @param raw_file_path complete file path to SNP data in traw format
#' @param window_list a list of window centers
#' @param window_size a numeric value indicating the size of windows in BP
#' @param pre_allocated_dir a directory where pre-allocated SNP window lists are kept
#'
#'
#' @return a list of lists, with each sub-list containing in its first element the SNP xposition and in the second element a Z matrix of SNPs (SKAT-appropriate format)
#' @export
#'
#' @examples
#' genotype_lists <- pre_allocate(raw_file_path = raw_file_path, window_list = window_list, window_size = window_size, window_shift = window_shift)
#'
# pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/SKAT_MC_parallel/pre_allocated"
pre_allocate <- function(raw_file_path, this_scaff_subset, window_size, window_shift,
                         pre_allocated_dir, window_list=NA){

  if (is.na(window_list)){
    message("Window list not provided to pre_allocate(), so we will create a list to run over the entire SNP data structure")
    if(!dir.exists(pre_allocated_dir)) dir.create(pre_allocated_dir,
                                                  recursive = TRUE)

    pre_allocated_path <- paste0(pre_allocated_dir, "/", basename(tools::file_path_sans_ext(raw_file_path)), ".", window_size, "bp_win.rds")

    if(file.exists(pre_allocated_path)){
      message(paste0("Data structure already exists. Read in from: ", pre_allocated_path))
      pos_and_SNP_list <- readRDS(pre_allocated_path)
    }
    if(!file.exists(pre_allocated_path)){
      this_scaff_subset <- data.table::fread(paste0(raw_file_path),
                                             fill = TRUE)

      message(paste0("Read genotype data for raw file ", raw_file_path,
                     " with dimensions", dim(this_scaff_subset)[1], "SNPs and ",
                     dim(this_scaff_subset)[2], "genotypes"))
      this_scaff_subset$POS <- as.numeric(as.character(this_scaff_subset$POS))

      max_window <- max(this_scaff_subset$POS)

      # The starting position is either A) one half of the window size,
      #   since we declare each window position to be the center of each window
      #   and wish to start at the start of the genome; or
      #   B) the first SNP position rounded to one-half of the window size,
      #   for situations in which the first window declared by approach A
      #   would not contain any SNP.
      starting_position <- max((window_size / 2),
                               plyr::round_any(
                                 min(this_scaff_subset$POS), window_size / 2))

      window_list <- seq(starting_position,
                         max_window,
                         by=window_shift)

      message("Allocating data structure")
      ptm <- proc.time()
      #browser()
      pos_and_SNP_list <- lapply(window_list,  function(x) extract_window(this_position = x,
                                                                          window_size = window_size,
                                                                          this_scaff_subset = this_scaff_subset))

      time <- proc.time() - ptm
      message(paste0("Took ", time[3], "s to allocate pos, SNP data structure"))

      saveRDS(pos_and_SNP_list, pre_allocated_path)
      }
  }

  # DEPRECATED? If so then remove this, and the prior if statement.
  # if (!is.na(window_list)){
  #   message("Allocating data structure")
  #   ptm <- proc.time()
  #   #browser()
  #   pos_and_SNP_list <- lapply(window_list,  function(x) extract_window(this_position = x,
  #                                                                       window_size = window_size,
  #                                                                       this_scaff_subset = this_scaff_subset))
  #
  #   time <- proc.time() - ptm
  #   message(paste0("Took ", time[3], "s to allocate pos, SNP data structure"))
  # }

  #browser()
  pos_and_SNP_list
}
