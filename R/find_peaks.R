#' Detect local maxima in time series
#'
#' @description
#' Detects peaks in a vector and calculates the peak height. This function is only appropriate for symmetric gaussian peaks and does not take into account any baseline correction as it required in 'real word' data. Therefore, it does not substitute sophisticated peak detection and integration tools and is only used for illustration purposes in our vignettes.
#'
#' @param df
#' A data frame containing x and y coordinates.
#'
#' @return A data frame containing x and y coordinates of peaks.
#'
#' @examples
#' ## create df
#' df <- data.frame(x = 1:1000, y = dnorm(1:1000,300,20))
#' ## plot
#' with(df, plot(x,y))
#' ## detect peak
#' find_peaks(df)
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
find_peaks <- function(df) {
    shape <- diff(sign(diff(df[["y"]], na.pad = FALSE)))
    pks <- which(shape < 0) # find potential turning points
    pks <- pks[which(df[["y"]][pks] > 0)] # remove zero intensity "peaks", they are chunk
    pks <- pks + 1
    if (any(diff(pks) == 1)) pks <- pks[-which(diff(pks) == 1)]
    return(df[pks,])
}#end find peaks
