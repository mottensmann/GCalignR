#' Simulating a simple chromatogram
#'
#'@param peaks vector of mean retention times representing peaks in a chromatogram
#'
#'@keywords internal
#'
simple_chroma <- function(peaks = c(10,13,25,37,50)) {
    x <- seq(from = 0, to = 60, length = 10000)
    y <- rep(0, length(x))
    for (i in 1:length(peaks)) y <- y + dnorm(x,mean = peaks[i], sd = sample(seq(0.2,0.75, 0.01),1))
    plot(x,y, type = "l", lwd = 1, xlab = "Retention time", ylab = "", axes = F)
    axis(side = 1,at = seq(0, 60, 10), labels = seq(0, 60, 10))
    return(data.frame(x = x, y = y))

}

