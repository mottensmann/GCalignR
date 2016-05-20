#' Determines the number of shared peaks between a chromatogram and a reference
#'
#' @description Estimate the most suitable linear shift in seconds to maximize the number of
#' shared peaks with the reference
#'
#' @inheritParams matrix_append
#' @inheritParams peak_shift
#'
#' @inheritParams align_chromatograms
#'
#' @param shift_steps
#' numeric vector of shifts for which the number of shared peaks is determined.
#'
#' @param error
#' numeric value defining the allowed difference in retention times in
#' of two shared peaks. The default \code{error=0} counts
#' a peak as shared if retention times match excatly.
#'
#'
#' @return
#' a data.frame containing the number of shared peaks and corresponding shifts for
#' all evaluated steps
#'
#' @details Shared peaks fall within the retention time of the reference and +- the Error (denoted in seconds)
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'
#' @export
shared_peaks <- function(gc_peak_df, ref_df, shift_steps, error=0, rt_col_name) {
    # Calculate the Number of shared peaks between a Chromatogram and its reference
    # Shared peaks fall within the retention time of the reference and +- the Error [s]
    ref <- ref_df[[rt_col_name]]
    no_of_peaks <- numeric(0)

    for (j in 1:length(shift_steps)){
        temp <- gc_peak_df[[rt_col_name]]+shift_steps[j] # Shift all Peaks by the same step
        peaks <- 0
        for (k in 1:length(temp)){ # loop through all Peaks of the current Sample
            for (l in 1:length(ref)){ # loop through the Reference Chromatogram
                temp_peak <- temp[k]
                ref_peak <- ref[l]
                if ( temp_peak!=0){ # Avoid comparison with cases of RT=0
                    if ((temp_peak <= ref_peak+error) & (temp_peak >= ref_peak-error)){
                        peaks <- peaks+1
                    }
                }
            }

        }
        no_of_peaks  <- c(no_of_peaks ,peaks)

    }
    output <- list(no_of_peaks ,shift_steps)
    return(output)
}
