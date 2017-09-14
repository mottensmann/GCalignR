#' GCalignR: A Package to Align Gas Chromatography Peaks Based on Retention Times
#'
#'@description
#' \strong{\code{GCalignR}} contains the functions listed below. Follow the links to access the documentation of each function.
#'
#'\code{\link{align_chromatograms}} executes all alignment steps.
#'
#'\code{\link{as.data.frame.GCalign}} exports aligned data to data frames.
#'
#'\code{\link{check_input}} tests the input data for formatting issues.
#'
#'\code{\link{draw_chromatogram}} visualises peak lists in form of a chromatogram.
#'
#'\code{\link{find_peaks}} detects and calculates peak heights in chromatograms. Not intended to be used for peak integration in empirical data. Used for illustration purposes only.
#'
#'\code{\link{gc_heatmap}} visualises aligned datasets using heatmaps that can be customised.
#'
#'\code{\link{norm_peaks}} allows to compute the relative abundance of peaks with samples.
#'
#'\code{\link{peak_interspace}} gives a histogram of the distance between peaks within samples over the whole dataset.
#'
#'\code{\link{read_peak_list}} reads the content of a text file and converts it to a list.
#'
#'\code{\link{remove_blanks}} removes peaks resembling contaminations from aligned datasets.
#'
#'\code{\link{remove_singletons}} removes peaks that are unique for one individual sample.
#'
#'\code{\link{simple_chroma}} creates simple chromatograms for testing and illustration purposes.
#'
#'@details
#' More details on the package are found in the vignettes that can be accessed via \code{browseVignettes("GCalignR")}.
#'
#' @docType package
#' @name GCalignR
#'
NULL
