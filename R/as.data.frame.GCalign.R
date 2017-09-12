#' Coearce aligned data to a data frame
#'
#' @description
#' Coverts aligned datasets into a data frame with columns for each peak and rows representing samples.
#'
#' @param x
#' An object of class "GCalign". See \code{\link{align_chromatograms}} for details.
#'
#' @inheritParams base::as.data.frame
#'
#' @return
#' A list of data frames for each variable of the dataset.
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
as.data.frame.GCalign <- function(x, row.names = NULL, optional = FALSE, ...) {
dat <- x[["aligned"]]
dat <- lapply(dat, function(y) y[-1])
dat <- lapply(dat, function(y) as.data.frame(t(do.call("cbind",y)), row.names, ...))
return(dat)
}
