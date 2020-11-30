#' Labeled window of alternative allele counts for Populus trichocarpa SNPs
#'
#' A dataset containing the position of a SNP window in the genome and all of
#' the alternative allele counts for each SNP within the SNP window
#'
#' Support for the Poplar GWAS dataset is provided by the U.S.
#' Department of Energy, Office of Science Biological and Environmental
#' Research (BER) via the Bioenergy Science Center (BESC) under Contract
#' No. DE-PS02-06ER64304.  The Poplar GWAS Project used resources of the
#' Oak Ridge Leadership Computing Facility and the Compute and Data Environment
#' for Science at Oak Ridge National Laboratory, which is supported by the
#' Office of Science of the U.S. Department of Energy under
#' Contract No. DE-AC05-00OR22725
#'
#' @format A list containing three values
#'
#' \describe{
#'   \item{Position}{Integer, the position of the center of a SNP window,
#'   in base pairs}
#'   \item{Z}{Matrix of integers containing alternative allele counts for
#'   each of 175 SNPs (columns) in the SNP window, across 882 genotypes (rows)}
#'   \item{Chr}{Integer or character indicating the chromosome or other
#'   scaffold, respectively, on which the SNP window is found}
#' }
#' @source \url{https://cbi.ornl.gov/data}
"sample_SNP_window"

#' Covariates for poplar samples in stem regeneration experiment
#'
#' A dataset containing covariates including stem diameter, phase and principal
#' components (PCs)
#'
#' @format A data frame with 882 rows and 11 variables
#'
#' \describe{
#'   \item{Stem_diameter}{Numeric value (double) representing the diameter of
#'   the stem of a poplar sample, in mm}
#'   \item{Phase1}{Indicator variable, with 1 if the sample was studied
#'   in Phase 1}
#'   \item{Phase2}{Indicator variable, with 1 if the sample was studied
#'   in Phase 2}
#'   \item{Phase3}{Indicator variable, with 1 if the sample was studied
#'   in Phase 3}
#'   \item{Phase4}{Indicator variable, with 1 if the sample was studied
#'   in Phase 4}
#'   \item{Phase5}{Indicator variable, with 1 if the sample was studied
#'   in Phase 5}
#'   \item{Phase6}{Indicator variable, with 1 if the sample was studied
#'   in Phase 6}
#'   \item{Phase7}{Indicator variable, with 1 if the sample was studied
#'   in Phase 7}
#'   \item{PC1}{Loadings for the first principal component derived from SNPs
#'   using the `pca` method in PLINK}
#'   \item{PC2}{Loadings for the second principal component derived from SNPs
#'   using the `pca` method in PLINK}
#'   \item{PC3}{Loadings for the third principal component derived from SNPs
#'   using the `pca` method in PLINK}
#'}
"sample_covariates"

#' SNP data in `.traw` format for poplar samples in stem regeneration experiment
#'
#' This dataset contains alternative allele counts for 5076 base positions from
#' about 14.46Mb 14.55Mb to on Chromosome 10 of Populus trichocarpa. This data
#' was obtained first by downloading
#' (\href{https://cbi.ornl.gov/data}{poplar SNP data}) from the Department
#' of Energy's Center for Bioenergy Innovation, converting the data from Variant
#' Call Format (VCF) to `.traw` using PLINK (see
#' \href{https://www.cog-genomics.org/plink2/formats}{PLINK documentation})
#' and finally loading the `.traw` file into R using
#' \code{\link[data.table]{fread}}. Each row represents a single nucleotide
#' polymorphism (SNP) with variation across the poplar GWAS population.
#'
#' Support for the Poplar GWAS dataset is provided by the U.S.
#' Department of Energy, Office of Science Biological and Environmental
#' Research (BER) via the Bioenergy Science Center (BESC) under Contract
#' No. DE-PS02-06ER64304.  The Poplar GWAS Project used resources of the
#' Oak Ridge Leadership Computing Facility and the Compute and Data Environment
#' for Science at Oak Ridge National Laboratory, which is supported by the
#' Office of Science of the U.S. Department of Energy under
#' Contract No. DE-AC05-00OR22725
#'
#' @format A data frame with 5076 rows and 888 variables
#'
#' \describe{
#'   \item{CHR}{Integer or character value indicating the chromosome or scaffold
#'   on which a given SNP is found}
#'   \item{SNP}{Character providing a label for a given SNP; this is optional
#'   and the column contains "." for all values by default if SNPs are not
#'   named. This variable is not used by the `mtmcskat` package.}
#'   \item{X.C.M}{Position of a given SNP in morgans or centimorgans; this is
#'   optional and can be filled with "0" if not used. This variable is not used
#'   by the `mtmcskat` package.}
#'   \item{POS}{Integer providing the position of a given SNP, in base pairs}
#'   \item{COUNTED}{Character, either "A", "T", "C" or "G" indicating the
#'   common (or most common) allele at the position of the given SNP}
#'   \item{COUNTER}{Character, either "A", "T", "C" or "G" indicating the
#'   alternative allele, also known as the rare allele,
#'   at the position of the given SNP. If multiple alternative alleles exist
#'   for a position, they are provided on separate rows with the same common
#'   allele and position.}
#'   \item{X201782_400194}{Columns 7-888 contain alternative allele
#'   counts for each of 882 genotypes in the poplar GWAS population. These are
#'   integer values ranging from 0-2, indicating the number of alleles which are
#'   the alternative allele for a given SNP in the column for a given genotype.
#'   The name of each column from 7-888 follows a format that includes the
#'   Family ID and Individual ID for each genotype, following the format
#'   X<FID>_<IID>}
#'}
"sample_genodata"

#' MTMCSKAT results from multithreading over SNP windows
#'
#' Results obtained from running MTMCSKAT with parallelization over SNP windows
#' in \code{\link{sample_re_allocated_SNP_windows}}, using the function
#' \code{\link{mtmcskat_SNPs}}
#'
#' @format A data frame with 3 rows and 4 variables
#' \describe{
#'   \item{Chr}{Integer or character value indicating the chromosome or scaffold
#'   on which a given SNP is found}
#'   \item{position}{Integer providing the position of a given SNP, in base
#'   pairs}
#'   \item{SKAT_p-val}{Numeric value (double) providing the p-value obtained by
#'   SKAT without resampling}
#'   \item{SKAT_p-val_resampled}{Numeric value (double) providing the empirical
#'   p-value obtained by MTMCSKAT after resampling}
#'}
"sample_mtmcskat_results"

#' MTMCSKAT results from multithreading over null models in a single batch
#'
#' Results obtained from running MTMCSKAT over
#' \code{\link{sample_re_allocated_SNP_windows}} with multithreading over
#' null models using the function \code{\link{mtmcskat_NullModels}},
#' while all null models can fit into the allotted RAM simultaneously and
#' thus all resampling can be completed in a single batch
#'
#' @format A data frame with 3 rows and 4 variables
#' \describe{
#'   \item{Chr}{Integer or character value indicating the chromosome or scaffold
#'   on which a given SNP is found}
#'   \item{position}{Integer providing the position of a given SNP, in base
#'   pairs}
#'   \item{SKAT_p-val}{Numeric value (double) providing the p-value obtained by
#'   SKAT without resampling}
#'   \item{SKAT_p-val_resampled}{Numeric value (double) providing the empirical
#'   p-value obtained by MTMCSKAT after resampling}
#'}
"sample_mtmcskatnm1_results"

#' MTMCSKAT results from multithreading over null models in two batches
#'
#' Results obtained from running MTMCSKAT over
#' \code{\link{sample_re_allocated_SNP_windows}} with multithreading over
#' null models using the function \code{\link{mtmcskat_NullModels}}, while not
#' all null models can fit into the allotted RAM simultaneously, thus requiring
#' resampling be completed in two batches
#'
#' @format A data frame with 3 rows and 4 variables
#' \describe{
#'   \item{Chr}{Integer or character value indicating the chromosome or scaffold
#'   on which a given SNP is found}
#'   \item{position}{Integer providing the position of a given SNP, in base
#'   pairs}
#'   \item{SKAT_p-val}{Numeric value (double) providing the p-value obtained by
#'   SKAT without resampling}
#'   \item{SKAT_p-val_resampled}{Numeric value (double) providing the empirical
#'   p-value obtained by MTMCSKAT after resampling}
#'}
"sample_mtmcskatnm2_results"

#' MTSKAT results from multithreading over SNP windows
#'
#' Results obtained from running MTSKAT with multithreading over SNP windows
#' in \code{\link{sample_re_allocated_SNP_windows}}, using the function
#' \code{\link{mtskat}}
#'
#' @format A data frame with 3 rows and 4 variables
#' \describe{
#'   \item{Chr}{Integer or character value indicating the chromosome or scaffold
#'   on which a given SNP is found}
#'   \item{position}{Integer providing the position of a given SNP, in base
#'   pairs}
#'   \item{SKAT_p-val}{Numeric value (double) providing the p-value obtained by
#'   SKAT without resampling}
#'   \item{SKAT_p-val_resampled}{Column containing NA, to be replaced with
#'   empirical p-values when this data structure is provided to
#'   \code{\link{mtmcskat_NullModels}} or \code{\link{mtmcskat_SNPs}}}
#'}
"sample_mtskat_results"

#' Tallies of p-values from resampling
#'
#' This dataset features tallies of permuted null model p-values that are above
#' or below the model p-value that was obtained from a null model without
#' resampling. The sums of p-values above or below this model p-value can be
#' used to calculate an empirical p-value. There is a row for each SNP window,
#' thread combination, since SNP windows are tested on all threads using
#' different null models on each thread. This sample dataset was obtained from
#' \code{\link{mtmcskat_NullModels}}, which multithreads Monte Carlo SKAT
#' over null models.
#'
#' @format A data frame with 192 rows and 4 variables
#' \describe{
#'   \item{position}{Integer indicating, in base pairs, the center of a SNP
#'   window}
#'   \item{SKAT_p-val}{Numeric value (double) indicating the original p-value
#'   for the given SNP window, obtained by SKAT without resampling}
#'   \item{n_perm_above}{The number of permuted p-values from a given thread
#'   that are above the original p-value for a given window}
#'   \item{n_perm_ltoreq}{The number of permuted p-values from a given thread
#'   that are below the original p-value for a given window}}
"sample_p_null_tallies"

#' Proportions of poplar image pixels labeled as shoot
#'
#' For all samples in the 882-genotype poplar GWAS population, this list
#' contains either 1) NA, or samples not studied in the sample GWAS, or 2)
#' a value ranging from 0 to 1, indicating the proportion of image pixels that
#' are labeled as shoot (rather than as callus or unregenerated stem). Samples
#' were treated by submerging poplar cuttings in water, applying TDZ to their
#' tips, and leaving a microcentrifuge tube ontop to provide humid conditions.
#' Phenotyping was performed five weeks after this treatment, via imaging plants
#' and classifying pixels with a trained convolutional neural network.
#'
#' @format A list with 882 values
"sample_phenotype"

#' List of consecutive SNP windows with positon labels
#'
#' This dataset contains a list of SNP windows, each formatted as described
#' in documentation for \code{\link{sample_SNP_window}}. All of these SNP
#' windows are consecutive and span the range of 14.46Mb 14.55Mb on Chr. 10.
#' Each SNP window spans a range of 3kb and consecutive windows are shifted by
#' 1kb relative to the previous adjacent window.
#'
#' @format A list with 90 values
"sample_pre_allocated_SNP_windows"

#' List of selected SNP windows significantly associated with shoot trait
#'
#' This dataset contains a list of SNP windows, each formatted as described
#' in documentation for \code{\link{sample_SNP_window}}. This list of SNP
#' windows includes all those within the range of 14.46Mb 14.55Mb on Chr. 10
#' that produced p-values between 0.01 and 0.001 in \code{\link{mtskat}}.
#' This list was produced by feeding
#' \code{\link{sample_pre_allocated_SNP_windows}} into
#' \code{\link{re_allocate_windows}} to yield a subset of SNP windows that
#' produced p-values within the desired range.
#'
#' @format A list with 4 values
"sample_re_allocated_SNP_windows"

#' List of SNP window centers
#'
#' This dataset is a list identical to that produced by the command
#' `seq(14460000, 14549000, by=1000)`. This is a list of SNP window positions,
#' where each value represents the center (in base pairs) of a SNP window.
#'
#' @format A list with 90 values
"sample_window_list"
