#' Full Alignment of Peak Lists by linear retention time correction.
#'
#' @description
#' Shifts all peaks within samples to maximise the similarity to a reference sample. For optimal results, a sufficient number of shared peaks are required to find a optimal solution. A reference needs to be specified, for instance using \code{\link{choose_optimal_reference}}. Linear shifts are evaluated within a user-defined window in discrete steps. The highest similarity score defines the shift that will be applied. If more than a single shift step yields to the same similarity score, the smallest absolute value wins in order to avoid overcompensation. The functions is envoked internally by \code{\link{align_chromatograms}}.
#'
#' @details
#' A similarity score is calculated as the sum of deviations in retention times between all reference peaks and the closest peak in the sample. The principle idea is that the appropriate linear transformation will reduce the deviation in retention time between homologous peaks, whereas all other peaks should deviate randomly. Among all considered shifts, the minimum deviation score is selected for subsequent full alignment by shifting all peaks of the sample by the same value.
#
#' @inheritParams align_chromatograms
#'
#' @inheritParams align_peaks
#'
#' @param reference
#' A character giving the name of a sample included in the dataset. All samples are aligned to the reference.
#'
#' @param step_size
#' Integer giving the step size in which linear shifts are evaluated between \code{max_linear_shift} and \code{-max_linear_shift}.
#'
#' @param Logbook
#' A list. If present, a summary of the applied linear shifts in full alignments of peak lists is appended to the list. If not specified, a list will be created automatically.
#'
#' @return
#' A list containing two items.
#' \item{chroma_aligned}{List containing the transformed data}
#' \item{Logbook}{Logbook, record of the applied shifts}
#'
#'
#' @examples
#' dat <- peak_data[1:10]
#' dat <- lapply(dat, function(x) x[1:50,])
#' x <- linear_transformation(gc_peak_list = dat, reference = "C2", rt_col_name = "time")
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
linear_transformation <- function(gc_peak_list,reference, max_linear_shift = 0.05, step_size = 0.01, rt_col_name, Logbook = NULL) {

    if (is.null(Logbook)) Logbook <- list()
    gc_peak_list <- remove_gaps(gc_peak_list, rt_col_name)
    # currently method is fixed
    method <- "Deviance"


    adjust_retention_time <- function(chromatogram, OptimalShift, ret_col_name){
        # Apply the best shift
        chromatogram[, ret_col_name] <- chromatogram[, ret_col_name] + OptimalShift
        return(chromatogram)
    }#end

    best_shift <- function(peaks, method) {
        # Determines the optimal shift
        shared <- as.vector(peaks[[1]])
        shifts <- as.vector(peaks[[2]])

        if (method == "Match") {
            # index of the best fit
            index <- which(shared == max(shared))
        } else if (method == "Deviance") {
            index <- which(shared == min(shared))
        }


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

    peak_shift <- function(gc_peak_df, ref_df, max_linear_shift = 0.05, step_size = 0.005, rt_col_name, method) {
        # 'peak_shift' uses 'shared_peaks' and 'best_shift' to find a suitable linear adjustment
        right_shift <- max_linear_shift
        left_shift <- max_linear_shift * -1
        shift_steps <- seq(from = left_shift ,to = right_shift,by = step_size)
        # Table of Shifts and shared Peaks
        output <- shared_peaks(gc_peak_df, ref_df, shift_steps, rt_col_name, method = method)
        # The Best shift
        output <- best_shift(peaks = output, method = method)
        return(output)
    }#end

    shared_peaks <- function(gc_peak_df, ref_df, shift_steps, rt_col_name, method)  {

        if (method == "Match") {
            #define function to count peaks
            num_sha <- function(ref, sam) { # for a given shift
                sum(unlist(lapply(X = ref, FUN = function(fx) {
                    ifelse(test = min(round(abs(fx - sam),2)) == 0, yes = 1, no = 0)
                })))
            }
            # apply
            no_of_peaks <- unlist(lapply(shift_steps, function(x) num_sha(ref = ref_df[[rt_col_name]], sam = gc_peak_df[[rt_col_name]] + x)))
                } else if (method == "Deviance") { #
                    # define function to estimate sumed differences
                    diff_sha <- function(ref, sam) { # for a given shift
                        sum(unlist(lapply(X = ref, FUN = function(fx) {
                            min(abs(fx - sam))
                        })))
                    }
                    # apply
                    no_of_peaks <- unlist(lapply(shift_steps, function(x) diff_sha(ref = ref_df[[rt_col_name]], sam = gc_peak_df[[rt_col_name]] + x)))

                }# end method
        output <- list(no_of_peaks, shift_steps)
        return(output)
    }#end

    shift_rts <- function(gc_peak_df, ref_df, max_linear_shift, step_size, method) {
        # Main Function doing the linear transformation.
        id <- gc_peak_df[["id"]][1]
        # drop the id column
        gc_peak_df <- gc_peak_df[-length(gc_peak_df)]
        optimal_shift <- peak_shift(gc_peak_df, ref_df, max_linear_shift, step_size, rt_col_name, method = method)
        shifted <- adjust_retention_time(gc_peak_df, optimal_shift, rt_col_name)
        # two lists per sample are created
        output <- list(shifted = shifted,optimal_shift = optimal_shift)
        return(output)
    }#end

    # ============================================================
if (max_linear_shift > 0) {
    ref <- gc_peak_list[[reference]]
    # temp list that allows to submit the name of the data.frames in lapply
    temp <- gc_peak_list
    id <- names(temp)
    for (j in 1:length(temp)) {
        temp[[j]][["id"]] <- id[j]
    }
    pbapply::pboptions(char = "+", style = 1)
    chrom_shift <- pbapply::pblapply(X = temp,FUN =  shift_rts, ref_df = ref, max_linear_shift = max_linear_shift, step_size = step_size, method = method)
    Logbook[["LinearShift"]] <- Logbooker(chrom_shift)
    chroma_aligned <- lapply(chrom_shift,function(x) x[-2])
} else {
    chroma_aligned <- gc_peak_list
    Logbook[["LinearShift"]] <- data.frame(shift = rep(0, length(gc_peak_list)), sample = names(gc_peak_list))
}
    return(list(chroma_aligned = chroma_aligned,Logbook = Logbook))
}#end linear_transformation

