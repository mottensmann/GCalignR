#' shifts peaks
#'
#' @description shifts peaks of individual peak data.frames.
#'
#' @param ref_df data.frame of form \code{gc_peak_df} serving as a reference.
#'
#' @inheritParams align_chromatograms
#'
#' @param step_size
#' defines the increments between \code{-max_linear_shift}:\code{max_linear_shift}
#' to consider in the search for an optimal linear transformation of peaks.
#'
#' @param error numeric value defining the allowed difference in retention times in
#'          derterming if two peaks are shared. The default \code{error=0} counts
#'          a peak a shared if retention times match excatly.
#'
#'
#' @return Numeric Value, indicating the best shift (e.g. -0.02 seconds)
#'
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'
#' @export


peak_shift <- function(gc_peak_df, ref_df, max_linear_shift=0.05, step_size=0.01, error=0, rt_col_name){
    # This functions shifts retention times of a chromatogram and estimates
    # the number of shared peaks with the reference.
    # Calling 'SharedPeaks' to count the number of shared peaks for a given shift
    # Calling 'BestShift' to find the best shift leading to the maximum similarity
    # If two Shifts lead to the maximu, take the adjustment with the smallest absolut value
    right_shift <- max_linear_shift
    left_shift <- max_linear_shift*-1
    shift_steps <- seq(from=left_shift ,to=right_shift,by=step_size)
    #     PeaksShared <- rep(0,length(ShiftSteps))
    #     PeaksLag <- rep(0,length(ShiftSteps))
    output <- shared_peaks(gc_peak_df, ref_df, shift_steps, error, rt_col_name) # List containg shared Peaks and their shifts
    output <- best_shift(output) # Which is the best setting
    output # Numeric Value, indicating the best shift (e.g. -0.02 seconds)
}
