% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blank_substraction.R
\name{blank_substraction}
\alias{blank_substraction}
\title{Subtraction of blank readings from sample readings}
\usage{
blank_substraction(input = NULL, blanks = NULL, conc_col_name = NULL)
}
\arguments{
\item{input}{A GCalign Object}

\item{blanks}{Character vector of names of negative controls.}

\item{conc_col_name}{If the input is a GCalign object the variable containing the abundance values needs to be specified.}
}
\description{
For each substance that is present in blanks, samples are corrected by subtraction of the respective quantity. If more than one sample is submitted, abundances are averaged. This procedure is sensitive to differences in the total concentration of samples and should be applied to samples where the preparation yields comparable concentrations for each sample.
}
\details{
Substances that are present in one or more blanks are identified in the aligned dataset, then the mean abundance is calculated for the blanks and the corresponding value is subtracted from each sample. If the control contains higher concentration (i.e. blank substation creates negative abundances) warnings will be shown and the respective value will be set to zero
}
\examples{
## Not run
#out <- blank_substraction(aligned_peak_data, blanks = "M2", conc_col_name = "area")


}
\keyword{beta}
