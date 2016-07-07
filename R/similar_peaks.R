#' Find peaks with similar retention time
#'
#' @description
#' Estimate degree of similarity between peaks by comparison of mean retention times of across samples.
#'
#' @param average_rts
#' vector of average retention times per peak across individuals.
#'
#' @inheritParams align_chromatograms
#'
#' @details
#' Similarity is evaluated at the level of min_diff_peak2peak, given in seconds.
#' The output is an interger vector containing indices of rows. Indices refer to a peak (i.e. substance)
#' that is similar to the previous peak For example index 5 indicates that the pair 4& 5 is similar.
#'
#' @return
#' integer vector of indices pointing to the last index of a row-pair that is similar
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @keywords internal
#'
similar_peaks <- function(average_rts, min_diff_peak2peak=0.05){

    difference <- rep(NA, (length(average_rts)-1)) # Estimate Difference between adjacent rows

    for (i in 2:length(average_rts)){
        difference[i] <- average_rts[i]-average_rts[i-1]
    }
    similar <- which(difference <= min_diff_peak2peak) # Which rows differ less the MinDistance ?
    return(similar)
}
