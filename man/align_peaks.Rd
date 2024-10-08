% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_peaks.R
\name{align_peaks}
\alias{align_peaks}
\title{align peaks individually among chromatograms}
\usage{
align_peaks(
  gc_peak_list,
  max_diff_peak2mean = 0.02,
  iterations = 1,
  rt_col_name,
  permute = TRUE,
  R = 1
)
}
\arguments{
\item{gc_peak_list}{List of data.frames. Each data.frame contains GC-data
(e.g. retention time, peak area, peak height) of one sample. Variables are
stored in columns. Rows represent distinct peaks. Retention time is a
required variable.}

\item{max_diff_peak2mean}{Numeric value defining the allowed deviation of the retention time of a given peak from the mean of the corresponding row (i.e. scored substance). This parameter reflects the retention time range in which peaks across samples are still matched as homologous peaks (i.e. substance). Peaks with retention times exceeding the threshold are sorted into a different row.}

\item{rt_col_name}{A character giving the name of the column containing the retention times. The decimal separator needs to be a point.}

\item{permute}{Boolean, by default a random permutation of samples is conducted prior for each row-wise alignment step. Setting this parameter to FALSE causes alignment of the dataset as it is.

order of samples is constantly randomised during the alignment. Allows to prevent this behaviour for maximal repeatability if needed.}

\item{R}{integer indicating the current iteration of the alignment step.
Created by \link{align_chromatograms}.}
}
\value{
a list of data.frames containing GC-data with aligned peaks.
}
\description{
\strong{align_peaks} allows to align similar peaks across samples
so that shared peaks are consistently located at the the same location (i.e.
defined as the same substance). The order of chromatograms (i.e. data.frames
in \code{gc_peak_list}) is randomized before each run of the alignment of
algorithm (if randomisation is not needed, this behaviour can be changed by setting \strong{permute = FALSE}). The main principle of this function is to reduce the variance in
retention times within rows, thereby peaks of similar retention time are
grouped together. Peaks that deviate significantly from the mean retention times
of the other samples are shifted to another row. At the start of a row the
first two samples are compared and separated if required, then all other
samples are included consecutively. If \code{iterations > 1} the whole
algorithm is repeated accordingly.
}
\details{
For each row the retention time of every sample is compared to the
mean retention time of all previously examined samples within the same row.
Starting with the second sample a comparison is done between the first and the
second sample, then between the third and the two first ones and so on.
Whenever the current sample shows a deviation from the mean retention time of
the previous samples a shift will either move this sample to the next row
(i.e. retention time above average) or all other samples will be moved to the
next row (i.e. retention time below average). If the retention time of the
sample in evaluation shows no deviation within \strong{-max_diff_peak2mean}:
\strong{max_diff_peak2mean} around the mean retention time no shifting is done
and the algorithm proceeds with the following sample.
}
\author{
Martin Stoffel (martin.adam.stoffel@gmail.com) & Meinolf Ottensmann
 (meinolf.ottensmann@web.de)
}
\keyword{internal}
