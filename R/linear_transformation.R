#' Shift peaks to eliminate systematic inaccuracies of peak detection by GC.
#'
#'@description
#'\code{linear_transformation()} applies small linear shifts of all peaks of
#'individual samples with respect to one reference. Thereby the number of shared
#'compounds among samples is maximized. The interval in which linear transformations are evaluated
#'is adjustable as well as the step size within this range.
#'
#' @param chromatograms list of data.frames with each being an individuals GC data.
#'
#' @param reference a character string indicating a sample included in \code{chromatograms}.
#'
#' @param shift the maximum shift in retention times that will be considered in estimating the best transformation.
#'
#' @param step_size indicates the step size in which linear shifts are evaluated
#'          between \code{shift} and \code{-shift}.
#'
#' @param error numeric value defining the allowed difference in retention times in
#'          derterming if two peaks are shared. The default \code{error=0} counts
#'          a peak a shared if retention times match excatly.
#'
#' @param rt_col_name character string denoting the column containing retention times
#'          in data.frames of \code{chromatograms}
#'
#' @return
#' A list of chromatograms with applied linear shifts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'

linear_transformation <- function(chromatograms,reference,
    shift=0.05, step_size=0.01, error=0, rt_col_name){
    # This is the master function which calls all sub-functions in order to
    # utilize a maximisation of the number of shared peaks
    # Mandatory arguments of this function are:
    # Chromatograms = List of Chromatograms, whereby each element of the List is a Matrix with the
    # peak extraction output (7 columns) of Xcalibur
    # References = Name(s) of Reference(s).

    # Include a vector of column names "ColNames" or specifiy the column which holds the
    # Apex of Retention Times "ColumnRT"

    ref <- chromatograms[[reference]]
    # Chroma_aligned <- list()

    shift_rts <- function(sample_df, ref_df, shift, step_size, error) {
        optimal_shift <- peak_shift(sample_df, ref_df, shift, step_size, error, rt_col_name)
        shifted <- adj_ret_time(sample_df, optimal_shift, rt_col_name)
    }

    chroma_aligned <- lapply(chromatograms, shift_rts, ref_df = ref, shift = shift, step_size = step_size, error = error)

    chroma_aligned
}

#' shifts peaks
#'
#' @description shifts peaks of individual chromatograms
#'
#' @param sample_df data.frame containing individual GC data
#'
#'  @param ref_df data.frame containing individual GC data of the reference chromatogram
#'
#' @param shift the maximum shift in retention times that will be considered in estimating the best transformation.
#'
#' @param step_size indicates the step size in which linear shifts are evaluated
#'          between \code{shift} and \code{-shift}.
#'
#' @param error numeric value defining the allowed difference in retention times in
#'          derterming if two peaks are shared. The default \code{error=0} counts
#'          a peak a shared if retention times match excatly.
#'
#' @param rt_col_name character string denoting the column containing retention times
#'          in data.frames of \code{chromatograms}
#'
#' @return Numeric Value, indicating the best shift (e.g. -0.02 seconds)
#'
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'
#' @export


peak_shift <- function(sample_df, ref_df, shift=0.05, step_size=0.01, error=0, rt_col_name){
    # This functions shifts retention times of a chromatogram and estimates
    # the number of shared peaks with the reference.
    # Calling 'SharedPeaks' to count the number of shared peaks for a given shift
    # Calling 'BestShift' to find the best shift leading to the maximum similarity
    # If two Shifts lead to the maximu, take the adjustment with the smallest absolut value
    right_shift <- shift
    left_shift <- shift*-1
    shift_steps <- seq(from=left_shift ,to=right_shift,by=step_size)
#     PeaksShared <- rep(0,length(ShiftSteps))
#     PeaksLag <- rep(0,length(ShiftSteps))
    output <- shared_peaks(sample_df, ref_df, shift_steps, error, rt_col_name) # List containg shared Peaks and their shifts
    output <- best_shift(output) # Which is the best setting
    output # Numeric Value, indicating the best shift (e.g. -0.02 seconds)
}


#' Estimate number of shared peaks between a chromatogram and a reference
#'
#' @description Estimate the most suitable linear shift in seconds to maximize the number of
#' shared peaks with the reference
#'
#' @param sample_df data.frame containing individual GC data
#'
#' @param ref_df data.frame containing individual GC data of the reference chromatogram
#'
#' @param shift_steps numeric vector containing steps to taken within the shift span
#'
#' @param error numeric value defining the allowed difference in retention times in
#'          derterming if two peaks are shared. The default \code{error=0} counts
#'          a peak a shared if retention times match excatly.
#'
#' @param rt_col_name character string denoting the column containing retention times
#'          in data.frames of \code{chromatograms}
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
shared_peaks <- function(sample_df, ref_df, shift_steps, error=0, rt_col_name) {
    # Calculate the Number of shared peaks between a Chromatogram and its reference
    # Shared peaks fall within the retention time of the reference and +- the Error [s]
    ref <- ref_df[[rt_col_name]]
    no_of_peaks <- numeric(0)

    for (j in 1:length(shift_steps)){
        temp <- sample_df[[rt_col_name]]+shift_steps[j] # Shift all Peaks by the same step
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

#' Find the best linear shift
#' @description Selects the optimal shifting time leading to the maximum number of shared peaks
#'
#' @param peak_list data.frame, the output from shared_peaks
#'
#' @return optimal shift in seconds to yield the maximum number of shared peaks
#'
#' @details If competing shifts exists (i.e. same number of peaks are shared),
#'          selects the smallest absolute value of a shift
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'
#' @export
#'
best_shift <- function(peak_list){

    ShiftTime <- as.vector(peak_list[[1]])
    Peaks <- as.vector(peak_list[[2]])
    Index <- which(ShiftTime==max(ShiftTime)) # find the best fit = maximum of shared peaks
    BestFit <- Peaks[Index]
    BestFit
    if (length(BestFit)>1){
        temp <- min(abs(BestFit)) # If equal shifts were found, take the smallest shift applied
        BestFit <- BestFit[BestFit==temp | BestFit==temp*-1]
        if (length(BestFit > 1)) BestFit <- BestFit[1]
    } else{
        BestFit
    }

}

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

adj_ret_time <- function(chromatogram, OptimalShift, ret_col_name){
    # Apply the estimated shift of the retention time to the
    # selected chromatogram to maximize the similarity compared to the reference
    chromatogram[, ret_col_name] <- chromatogram[, ret_col_name] + OptimalShift
    chromatogram
}
