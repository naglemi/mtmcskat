% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow_utils.R
\name{determine_subset_pval_boundaries}
\alias{determine_subset_pval_boundaries}
\title{Determine boundaries for subsetting SNP windows based on intial p-values}
\usage{
determine_subset_pval_boundaries(
  leading_0s,
  desired_sig_figs,
  terminal_resampling
)
}
\arguments{
\item{leading_0s}{An integer for the number of 0's between the decimal point
and the first nonzero value in the maximum p-value for SNP windows that will
be included in the subset of interest}

\item{desired_sig_figs}{An integer for the number of significant figures
the user desires for empirical p-values}

\item{terminal_resampling}{If `TRUE`, the lower boundary will be 0 so all
remaining SNP windows will be included in the final set for resampling, even
if their empirical p-values cannot be calculated out to the user-specified
number of significant figures}
}
\value{
A list with two decimal values, named `upper` and `lower`
}
\description{
Determine boundaries for subsetting SNP windows based on intial p-values
}
\examples{

# Produce boundaries from 0.01 to 0.001
determine_subset_pval_boundaries(leading_0s = 2,
                                 desired_sig_figs = 2,
                                 terminal_resampling = FALSE)

# Produce boundaries from 0.01 to 0
determine_subset_pval_boundaries(leading_0s = 2,
                                 desired_sig_figs = 2,
                                 terminal_resampling = TRUE)

# Produce boundaries from 0.01 to 0.0001
determine_subset_pval_boundaries(leading_0s = 2,
                                 desired_sig_figs = 3,
                                 terminal_resampling = FALSE)

}
