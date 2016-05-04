#' maximisation of the number of shared peaks with reference chromatogram
#' 
#' @param chromatograms list of data.frames with each data.frame being an individuals gc data
#' @param high Upper threshold for retention times. RTs lower than \code{high} will be kept
#' @param rt_col_name character string for name of retention time in gc table
#' 
#' @return 
#' gc table with retention times cut
#'
#' @references 
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
        shifted <- adj_ret_time(sample_df, optimal_shift) 
    }
    
    chroma_aligned <- lapply(chromatograms, shift_rts, ref_df = ref, shift = shift, step_size = step_size, error = error)
    
    chroma_aligned
}

#' shifts retention times of a chromatogram and estimates the number of shared peaks with the
#' reference.
#' 
#' @param sample_df \code{data.frame} with individual gc data
#' @param ref_df \code{data.frame} with individual gc data of the reference chromatogram
#' @param shift shift span to try for maximising shared peaks
#' @param step_size steps to take within the shift span
#' @param error allowed error
#' @param rt_col_name RT time column name
#'   
#' @return Numeric Value, indicating the best shift (e.g. -0.02 seconds)
#' 
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'   
 

#### doc missing
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


#' Calculate the Number of shared peaks between a Chromatogram and its reference
#' 
#' 
#' @param sample_df \code{data.frame} with individual gc data
#' @param ref_df \code{data.frame} with individual gc data of the reference chromatogram
#' @param shift_steps steps to taken within the shift span
#' @param error allowed error
#' @param rt_col_name RT time column name
#'   
#' @return list of number of peaks and shift steps
#' 
#' @details Shared peaks fall within the retention time of the reference and +- the Error [s]
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'   
#'   
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
}

#' Selects the optimal Shifting time leading to the maximum number of shared peaks
#' 
#' @param list, output from shared_peaks
#'   
#' @return list of number of peaks and shift step for which this is maximised
#' 
#' @details If competing shifts exists (i.e. same number of peaks are shared), 
#'          selects the smallest absolute value of a shift
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'   
#'  
#'  
best_shift <- function(list){

    ShiftTime <- as.vector(list[[1]])
    Peaks <- as.vector(list[[2]])
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

#' Apply the estimated shift of the retention time to the selected chromatogram 
#' to maximize the similarity compared to the reference
#' 
#' @param list, output from shared_peaks
#'   
#' @return list of number of peaks and shift step for which this is maximised
#' 
#' @details If competing shifts exists (i.e. same number of peaks are shared), 
#'          selects the smallest absolute value of a shift
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'   
#'  

adj_ret_time <- function(Data, OptimalShift, ret_col_name){
    # Apply the estimated shift of the retention time to the 
    # selected chromatogram to maximize the similarity compared to the reference
    Data[, ret_col_name] <- Data[, ret_col_name] + OptimalShift
    Data
}
