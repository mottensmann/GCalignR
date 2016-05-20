#'Check whether similar rows are redundant or not
#'
#' @description
#' \code{check_redundancy()} evaluates if peaks of similar retention time are redundant.
#' This is ture as long as any sample contains peaks in either one or none of the
#' row-pair, but never in both.
#'
#' @inheritParams matrix_append
#'
#' @inheritParams align_chromatograms
#'
#' @param similar_peaks
#' integer indexing the position of a peaks that is similar to the one before. Output from \code{similar_peaks}.
#'
#'
#' @return
#' a binary vector of ones and zeros. One indicates redundancy.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'

check_redundancy <- function(gc_peak_df, similar_peaks, rt_col_name){
    # If only one of two neighbouring rows contain a substance
    # they are redundant, coded by a One
    row1 <- gc_peak_df[similar_peaks-1, rt_col_name] # Extract previous row
    row2 <- gc_peak_df[similar_peaks, rt_col_name] # Extract current row
    redundant <- 0
    if (row1==0 | row2==0){
        redundant <- 1
    }
    return(redundant)
}
