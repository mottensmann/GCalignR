#' Shift peaks to eliminate systematic inaccuracies of peak detection by GC.
#'
#'@description
#'\code{linear_transformation()} applies small linear shifts of all peaks of
#'individual samples with respect to one reference. Thereby the number of shared
#'compounds among samples is maximized. The interval in which linear transformations are evaluated
#'is adjustable as well as the step size within this range.
#'
#'
#' @param reference
#' a character string indicating a sample included in \code{gc_peak_list} used as a reference to align to.
#'
#' @param step_size
#' indicates the step size in which linear shifts are evaluated
#' between \code{max_linear_shift} and \code{-max_linear_shift}.
#'
#' @inheritParams align_chromatograms
#'
#' @inheritParams align_peaks
#'
#' @return
#' A list of data.frames containing chromatograms with applied linear shifts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @keywords internal
#' @export
#'

linear_transformation <- function(gc_peak_list,reference,
                                  max_linear_shift=0.05, step_size=0.01, rt_col_name,Logbook){

    # Defining internal functions #
    ###############################
    # peak_shift() determines the optimal shift by use of shared_peaks() and best_shift()
    peak_shift <- function(gc_peak_df, ref_df, max_linear_shift=0.05, step_size=0.005, error=0, rt_col_name){
        right_shift <- max_linear_shift
        left_shift <- max_linear_shift*-1
        shift_steps <- seq(from=left_shift ,to=right_shift,by=step_size)
        output <- shared_peaks(gc_peak_df, ref_df, shift_steps, error, rt_col_name) # Table of Shifts and shared Peaks
        output <- best_shift(output) # Best shift
        return(output)
    }

    shared_peaks <- function(gc_peak_df, ref_df, shift_steps, error=0, rt_col_name) {
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

    best_shift <- function(peaks){
# Which is the appropriate shift to select?
        shared <- as.vector(peaks[[1]])
        shifts <- as.vector(peaks[[2]])
        index <- which(shared==max(shared)) # index of the best fit
        BestFit <- shifts[index] # Best Value

        if (length(BestFit)>1){
            temp <- min(abs(BestFit)) # If equal shifts were found, take the smallest shift applied
            BestFit <- BestFit[BestFit==temp | BestFit==temp*-1]
            if (length(BestFit > 1)) BestFit <- BestFit[1]
        }else{
            # BestFit
        }
#         if(file.exists(paste0(strsplit(as.character(match.call(definition = sys.function(sys.parent(5)), call = sys.call(sys.parent(5)))["data"]),split = ".txt"),"_LogFile.txt"))){
#             sink(paste0(strsplit(as.character(match.call(definition = sys.function(sys.parent(5)), call = sys.call(sys.parent(5)))["data"]),split = ".txt"),"_LogFile.txt"),append = TRUE)
#             cat(paste('\nShift = ',as.character(format(round(BestFit,3),nsmall=2)),'\tShared Peaks = ',as.character(shared[index[1]]))) # Delete later
#             sink()
#         }
        return(BestFit)
    }
    adjust_retention_time <- function(chromatogram, OptimalShift, ret_col_name){
        # Apply the best shift
        chromatogram[, ret_col_name] <- chromatogram[, ret_col_name] + OptimalShift
        return(chromatogram)
    }

    ###############################
#
#     if(file.exists(paste0(strsplit(as.character(match.call(definition = sys.function(sys.parent(1)), call = sys.call(sys.parent(1)))["data"]),split=".txt"),"_LogFile.txt"))){
#         sink(paste0(strsplit(as.character(match.call(definition = sys.function(sys.parent(1)), call = sys.call(sys.parent(1)))["data"]),split=".txt"),"_LogFile.txt"),append = TRUE)
#         cat("\nSamples in order of comparisons with the reference:\n")
#         print(names(gc_peak_list))
#         sink()
#     }

    ref <- gc_peak_list[[reference]]

    shift_rts <- function(gc_peak_df, ref_df, max_linear_shift, step_size, error) {
    # Main Function doing the linear transformation.
    id <- gc_peak_df[["id"]][1]
    gc_peak_df <- gc_peak_df[-length(gc_peak_df)] # drop the id column
    optimal_shift <- peak_shift(gc_peak_df, ref_df, max_linear_shift, step_size, error, rt_col_name)
    shifted <- adjust_retention_time(gc_peak_df, optimal_shift, rt_col_name)
    # Logbook[["LinearShift"]][id] <- as.character(optimal_shift)
    output <- list(shifted=shifted,optimal_shift=optimal_shift) # two lists per sample are created
    return(output)
    }

    temp <- gc_peak_list # temp list that allow to submit the name of the data.frames in lapply
    id <- names(temp)
    for(j in 1:length(temp)){
        temp[[j]][["id"]] <- id[j]
    }
    Logbooker <- function(chrom_shift){# to extract optimal shifts applied
        shift <- unlist(lapply(chrom_shift, function(x) x[-1]))
        sample <- strsplit(names(shift),split = ".optimal_shift")
        sample <- unlist(sample)
        shift <- as.vector(shift)
        out <- data.frame(shift,sample)
        return(out)
    }
    chrom_shift <- lapply(X = temp,FUN =  shift_rts, ref_df = ref, max_linear_shift = max_linear_shift, step_size = step_size, error = 0)
    Logbook[["LinearShift"]] <- Logbooker(chrom_shift)
    chroma_aligned <- lapply(chrom_shift,function(x) x[-2])
    return(list(chroma_aligned=chroma_aligned,Logbook=Logbook))
}#end linear transform

correct_colnames <- function(gc_peak_df,col_names){
    colnames(gc_peak_df) <- col_names
    return(gc_peak_df)
}
