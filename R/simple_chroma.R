#' Simulating a simple chromatogram
#'
#'@param peaks vector of mean retention times representing peaks in a chromatogram
#'
#'@keywords internal
#'
simple_chroma <- function(peaks = c(10,13,25,37,50), N = 1, min = 0, max = 30) {


    ## internal functions
    ## ##################

    lin_error <- function(range = 5) sample(x = seq(from = range*-1, to = range, by = 0.01),size = 1)
    rand_error <- function(range = 0.05, peaks = NULL) sample(x = seq(from = range*-1, to = range, by = 0.01), size = length(peaks),replace = T)

    # function creating single chromatograms
    fx <- function(lin_size = NULL) {
        if (N > 1) {
        peaks <- sample(x = peaks, size = sample(x = round((length(peaks)*0.8)):length(peaks), size = 1))
        peaks <- peaks + lin_size
        peaks <- peaks + sample(x = seq(from = -0.02, by = 0.01, to = 0.02), size = length(peaks), replace = T, prob = c(0.05,0.15,0.6,0.15,0.05))
        }

        for (i in 1:length(peaks)) y <- y + dnorm(x,mean = peaks[i], sd = sample(seq(0.2,0.4, 0.01),1))
        return(y)
    }
    ## #################

    #if (length(peaks) < 10) stop("Minimum population size is 3 peaks")
    # vector of retention times
    x <- seq(from = min, to = max, length = 10000)
    # vector for intensities
    y <- rep(0, length(x))

    # data frame to store simulated chromatogram data
    df <- data.frame(x = rep(x, N), y = rep(y, N), sample = rep(paste0("A", as.character(1:N)), each = length(x)))

    # preallocate a vector to write inensities to
    y2 <- numeric(0)



    # for all N´s
    for (i in 1:N) y2 <- c(y2, fx(lin_size = i*0.7))

    # update data frame
    df[["y"]] <- y2

    return(df)
}

find_peaks <- function(df) {
    shape <- diff(sign(diff(df[["y"]], na.pad = FALSE)))
    pks <- which(shape < 0) # find potential turning points
    pks <- pks[which(df[["y"]][pks] > 0)] # remove zero intensity "peaks", they are chunk
    pks <- pks + 1
    if (any(diff(pks) == 1)) pks <- pks[-which(diff(pks) == 1)]
    return(df[pks,])
}


