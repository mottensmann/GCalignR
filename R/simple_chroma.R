#' Simulate simple chromatograms
#'
#' @description
#' Creates chromatograms with user defined peaks for illustrative purposes. Linear drift is applied in sample order if more than one sample is created. See parameters of the function.
#'
#' @param peaks
#' A numeric vector giving the retention times on which gaussian distribution, defining peaks, are centered. If more than one sample is generated \code{N > 1}, \code{peaks} defines a population of peaks, from which samples are generated.
#'
#' @param N
#' An integer giving the number of chromatograms to create. By default \code{N = 1}.
#'
#' @param min
#' A numeric giving the minimum retention time.
#'
#' @param max
#' A numeric giving the maximum retention time.
#'
#' @param Names
#' A character vector giving sample names. If not specified, names are generated automatically.
#'
#' @param sd
#' A numeric vector of the same length as peaks giving the standard deviation of each peak. Only supported if N = 1.
#'
#' @return A data frame containing x and y coordinates and corresponding sample names.
#'
#' @examples
#' ## create a chromatogram
#' x <- simple_chroma(peaks = c(5,10,15), N = 1, min = 0, max = 30, Names = "MyChroma")
#' ## plot chromatogram
#' with(x, plot(x,y, xlab = "time", ylab = "intensity"))
#'
#'@author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
simple_chroma <- function(peaks = c(10,13,25,37,50), N = 1, min = 0, max = 30, Names = NULL, sd = NULL) {
if (is.null(Names)) Names <- paste0("A", as.character(1:N))
if (length(Names) != N) stop("Length of Names != N")

    ## internal functions
    ## ##################

    # lin_error <- function(range = 5) sample(x = seq(from = range*-1, to = range, by = 0.01),size = 1)
    # rand_error <- function(range = 0.05, peaks = NULL) sample(x = seq(from = range*-1, to = range, by = 0.01), size = length(peaks),replace = T)

    # function creating single chromatograms
    fx <- function(lin_size = NULL) {
        if (N > 1) {
        peaks <- sample(x = peaks, size = sample(x = round((length(peaks)*0.8)):length(peaks), size = 1))
        peaks <- peaks + lin_size
        peaks <- peaks + sample(x = c(-0.4,-0.2,0,0.2,0.4), size = length(peaks), replace = T, prob = c(0.05,0.15,0.6,0.15,0.05))
        }
if (N == 1 & !is.null(sd)) {
    std <- sd
    for (i in 1:length(peaks)) y <- y + dnorm(x,mean = peaks[i], sd = std[i])
} else {
    for (i in 1:length(peaks)) y <- y + dnorm(x,mean = peaks[i], sd = sample(seq(0.2,0.4, 0.01),1))
}
        return(y)
}#end fx

    # vector of retention times
    x <- seq(from = min, to = max, length = 10000)
    # vector for intensities
    y <- rep(0, length(x))

    # data frame to store simulated chromatogram data
    df <- data.frame(x = rep(x, N), y = rep(y, N), sample = rep(Names, each = length(x)))

    # preallocate a vector to write inensities to
    y2 <- numeric(0)

    # for all Ns
    for (i in 1:N) y2 <- c(y2, fx(lin_size = i*0.7))

    # update data frame
    df[["y"]] <- y2

    return(df)
}#end simple chroma



