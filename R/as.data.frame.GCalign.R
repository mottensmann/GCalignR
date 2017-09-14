#' Output aligned data in form of a data frame for each variable
#'
#' @description
#' Based on an object of class "GCalign" that was created using \code{\link{align_chromatograms}}, a list of data frames for each variable in the dataset is returned. Within data frames rows represent substances and columns are variables (i.e. substances).
#'
#' @param x
#' An object of class "GCalign". See \code{\link{align_chromatograms}} for details.
#'
#' @inheritParams base::as.data.frame
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @examples
#' data("aligned_peak_data")
#' out <- as.data.frame(x = aligned_peak_data)
#'
#' @export
#'
as.data.frame.GCalign <- function(x, row.names = NULL, optional = FALSE, ...) {
dat <- x[["aligned"]]
dat <- lapply(dat, function(y) y[-1])
dat <- lapply(dat, function(y) as.data.frame(t(do.call("cbind",y)), row.names, ...))
return(dat)
}
