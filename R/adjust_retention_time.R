#' Apply a linear shift to all retention times of a chromatogram
#' @description Shifts retention times of a sample by an previously defined value to optimise the
#'              similarity to the reference list by maximizing the number of shared  peaks.
#'
#' @param chromatogram data.frame containing gc data of an individual sample.
#'
#' @param OptimalShift numeric value indicating the optimal shift to apply to
#'          maximize the similarity to a reference.
#' @param ret_col_name character string denoting the column containing retention times in \code{chromatogram}.
#'
#' @return
#' \code{chromatogram}{data.frame cotaining linear adjusted chromatogram}
#'
#' @details If competing shifts exists (i.e. same number of peaks are shared),
#'          selects the smallest absolute value of a shift.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'
#' @export

adjust_retention_time <- function(chromatogram, OptimalShift, ret_col_name){
    # Apply the estimated shift of the retention time to the
    # selected chromatogram to maximize the similarity compared to the reference
    chromatogram[, ret_col_name] <- chromatogram[, ret_col_name] + OptimalShift
    chromatogram
}
