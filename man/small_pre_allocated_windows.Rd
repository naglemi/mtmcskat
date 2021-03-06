% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{small_pre_allocated_windows}
\alias{small_pre_allocated_windows}
\title{List of consecutive SNP windows with positon labels}
\format{
A list with 90 values
}
\usage{
small_pre_allocated_windows
}
\description{
This dataset contains a list of SNP windows, each formatted as described
in documentation for \code{\link{sample_SNP_window}}. All of these SNP
windows are consecutive and span the range of 14.49Mb 14.52Mb on Chr. 10.
Each SNP window spans a range of 3kb and consecutive windows are shifted by
1kb relative to the previous adjacent window.
}
\details{
This object has not undergone processing to filter out SNPs with high
missing rates or to impute missing values to the mean. Thus, the use of this
object will produce several warning messages from `SKAT`, noting that
ad-hoc missing rate filtering and imputation is being performed.
}
\keyword{datasets}
