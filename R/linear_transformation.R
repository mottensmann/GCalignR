#' Shift peaks to eliminate systematic inaccuracies of peak detection by GC.
#'
#'@description
#'\strong{linear_transformation()} shifts all peaks within chromatograms to maximise the number
#'of shared peaks with a reference chromatogram. Optimally, the reference contains known peaks which
#'also occur in the samples. If a sample is taken as a reference, samples with high concentrations
#'and clear peaks will lead to a better result.
#'
#' @param reference
#' character string with the name of a sample included in \code{gc_peak_list} used as a reference to align to.
#'
#' @param step_size
#' integer, indicating the step size in which linear shifts are evaluated between \strong{max_linear_shift} and \strong{-max_linear_shift}.
#'
#' @inheritParams align_chromatograms
#'
#' @inheritParams align_peaks
#'
#' @return
#' \item{chroma_aligned}{Transformed data}
#' \item{Logbook}{Logbook, record of the applied shifts}
#' List of data.frames containing chromatograms with applied linear shifts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @keywords internal
#'
linear_transformation <- function(gc_peak_list,reference,max_linear_shift = 0.05, step_size = 0.01, rt_col_name,Logbook){

# Revision 24-04-2017
# Prior to this date, all values were rounded to two decimals within the main function align_chromatograms. Now, this step is subsituted by  calculations based on rounded values within this function.
# Reordered alphabetically

    # Defining internal functions, The main function is peak_shift
    # ============================================================

    adjust_retention_time <- function(chromatogram, OptimalShift, ret_col_name){
        # Apply the best shift
        chromatogram[, ret_col_name] <- chromatogram[, ret_col_name] + OptimalShift
        return(chromatogram)
    }#end

    best_shift <- function(peaks) {
        # Determines the optimal shift
        shared <- as.vector(peaks[[1]])
        shifts <- as.vector(peaks[[2]])
        # index of the best fit
        index <- which(shared == max(shared))
        # Best Value
        BestFit <- shifts[index]
        if (length(BestFit) > 1) {
            # If equal shifts were found, take the smallest shift applied
            temp <- min(abs(BestFit))
            BestFit <- BestFit[BestFit == temp | BestFit == temp * -1]
            if (length(BestFit > 1)) BestFit <- BestFit[1]
        }else{
            # BestFit
        }
        return(BestFit)
    }#end

    correct_colnames <- function(gc_peak_df,col_names) {
        colnames(gc_peak_df) <- col_names
        return(gc_peak_df)
    }#end

    # to extract optimal shifts applied
    Logbooker <- function(chrom_shift) {
        shift <- unlist(lapply(chrom_shift, function(x) x[-1]))
        sample <- strsplit(names(shift),split = ".optimal_shift")
        sample <- unlist(sample)
        shift <- as.vector(shift)
        out <- data.frame(shift,sample)
        return(out)
    }#end

    peak_shift <- function(gc_peak_df, ref_df, max_linear_shift = 0.05, step_size = 0.005, error = 0, rt_col_name) {
        # 'peak_shift' uses 'shared_peaks' and 'best_shift' to find a suitable linear adjustment
        right_shift <- max_linear_shift
        left_shift <- max_linear_shift * -1
        shift_steps <- seq(from = left_shift ,to = right_shift,by = step_size)
        # Table of Shifts and shared Peaks
        output <- shared_peaks(gc_peak_df, ref_df, shift_steps, error, rt_col_name)
        # The Best shift
        output <- best_shift(output)
        return(output)
    }#end

    shared_peaks <- function(gc_peak_df, ref_df, shift_steps, error = 0, rt_col_name) {
        # get the rts of the reference
        ref <- ref_df[[rt_col_name]]
        # preallocate a vector to estimate the number of shared peaks
        no_of_peaks <- numeric(0)
        for (j in 1:length(shift_steps)) {
            # Shift all Peaks by the same step
            temp <- gc_peak_df[[rt_col_name]] + shift_steps[j]
            peaks <- 0
            # loop through all Peaks of the current Sample
            for (k in 1:length(temp)) {
                # loop through the Reference Chromatogram
                for (l in 1:length(ref)) {
                    temp_peak <- temp[k]
                    ref_peak <- ref[l]
                    # Avoid comparison with cases of RT = 0
                    if ( temp_peak != 0) {
                        if ((round(temp_peak, 2) <= round(ref_peak, 2) + error) & (round(temp_peak, 2) >= round(ref_peak, 2) - error)) {
                            peaks <- peaks + 1
                        }
                        # if ((temp_peak <= ref_peak + error) & (temp_peak >= ref_peak - error)) {
                        #     peaks <- peaks + 1
                        # }
                    }
                }
            }
            no_of_peaks  <- c(no_of_peaks, peaks)
        }
        output <- list(no_of_peaks, shift_steps)
        return(output)
    }#end

    shift_rts <- function(gc_peak_df, ref_df, max_linear_shift, step_size, error) {
        # Main Function doing the linear transformation.
        id <- gc_peak_df[["id"]][1]
        # drop the id column
        gc_peak_df <- gc_peak_df[-length(gc_peak_df)]
        optimal_shift <- peak_shift(gc_peak_df, ref_df, max_linear_shift, step_size, error, rt_col_name)
        shifted <- adjust_retention_time(gc_peak_df, optimal_shift, rt_col_name)
        # two lists per sample are created
        output <- list(shifted = shifted,optimal_shift = optimal_shift)
        return(output)
    }#end

    # ============================================================

    ref <- gc_peak_list[[reference]]
    # temp list that allows to submit the name of the data.frames in lapply
    temp <- gc_peak_list
    id <- names(temp)
    for (j in 1:length(temp)) {
        temp[[j]][["id"]] <- id[j]
    }
    chrom_shift <- lapply(X = temp,FUN =  shift_rts, ref_df = ref, max_linear_shift = max_linear_shift, step_size = step_size, error = 0)
    Logbook[["LinearShift"]] <- Logbooker(chrom_shift)
    chroma_aligned <- lapply(chrom_shift,function(x) x[-2])
    return(list(chroma_aligned = chroma_aligned,Logbook = Logbook))
}#end linear_transformation

