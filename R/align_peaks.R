#'align peaks individually among chromatograms
#'
#'@description \strong{align_peaks} allows to align similar peaks across samples
#'so that shared peaks are consistently located at the the same location (i.e.
#'defined as the same substance). The order of chromatograms (i.e. data.frames
#'in \code{gc_peak_list}) is randomized before each run of the alignment of
#'algorithm (if randomisation is not needed, this behaviour can be changed by setting \strong{permute = FALSE}). The main principle of this function is to reduce the variance in
#'retention times within rows, thereby peaks of similar retention time are
#'grouped together. Peaks that deviate significantly from the mean retention times
#'of the other samples are shifted to another row. At the start of a row the
#'first two samples are compared and separated if required, then all other
#'samples are included consecutively. If \code{iterations > 1} the whole
#'algorithm is repeated accordingly.
#'
#'@details For each row the retention time of every sample is compared to the
#'mean retention time of all previously examined samples within the same row.
#'Starting with the second sample a comparison is done between the first and the
#'second sample, then between the third and the two first ones and so on.
#'Whenever the current sample shows a deviation from the mean retention time of
#'the previous samples a shift will either move this sample to the next row
#'(i.e. retention time above average) or all other samples will be moved to the
#'next row (i.e. retention time below average). If the retention time of the
#'sample in evaluation shows no deviation within \strong{-max_diff_peak2mean}:
#'\strong{max_diff_peak2mean} around the mean retention time no shifting is done
#'and the algorithm proceeds with the following sample.
#'
#'@param gc_peak_list List of data.frames. Each data.frame contains GC-data
#'  (e.g. retention time, peak area, peak height) of one sample. Variables are
#'  stored in columns. Rows represent distinct peaks. Retention time is a
#'  required variable.
#'
#'@inheritParams align_chromatograms
#'
#'@param R integer indicating the current iteration of the alignment step.
#'  Created by \link{align_chromatograms}.
#'
#'@return a list of data.frames containing GC-data with aligned peaks.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@keywords internal
#'
align_peaks <- function(gc_peak_list, max_diff_peak2mean = 0.02, iterations = 1, rt_col_name, permute = TRUE, R = 1) {
    # print to Console
    # cat(paste('Iteration',as.character(R),'out of',as.character(iterations),' ... '))

    # Define the start
    current_row <- 1
    while (current_row != 'stop') {

        ##############
        ### timer ####
        ##############
        total <- nrow(gc_peak_list[[1]])
        # create progress bar

        if (interactive()) {
        pb <- utils::txtProgressBar(min = 0, max = total, style = 3, char = "+", width = 80)
        utils::setTxtProgressBar(pb, ifelse(is.numeric(current_row), current_row, total))
        }
        ##############


        # Randomize the order of samples
        if (isTRUE(permute)) {
            shuffle_order <- sample(1:length(gc_peak_list))
            gc_peak_list <- gc_peak_list[shuffle_order]
        }

        # Start with the second sample (S = 2)
    for (s in 2:length(gc_peak_list)) {
        # Select the rt of s
        current_rt <- gc_peak_list[[s]][current_row, rt_col_name]
        av_rt <- mean_retention_time_row(gc_peak_list, samples = 1:(s - 1), retention_row = current_row, rt_col_name)
        # Do nothing, if the selected is the only peak of that row
    if (is.na(av_rt)) {
        av_rt <- current_rt
    }
        # If the current sample is empty in that row
    if (current_rt == 0) {
        current_rt <- av_rt
    }
        # Current rt is larger than the mean
    if (current_rt > av_rt + max_diff_peak2mean) {
        gc_peak_list <- shift_rows(gc_peak_list,s,current_row)
        # current rt is smaller than the mean
    } else if (current_rt < (av_rt - max_diff_peak2mean)) {
        for (j in 1:(s - 1)) {
            if (gc_peak_list[[ j ]][current_row, rt_col_name] <= (current_rt + max_diff_peak2mean)) {
              # Do nothing, substance? position is valid
            } else {
                gc_peak_list <- shift_rows(gc_peak_list, j ,current_row)
            }
        }
    }
    }

    # Got to the next row
    current_row <- current_row + 1
    # Equalise the matrix size
    gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)
    last_substance_index <- length(mean_retention_times(gc_peak_list, rt_col_name)[mean_retention_times(gc_peak_list, rt_col_name) > 0])
    # remove unused rows
    gc_peak_list <- lapply(gc_peak_list, function(x) x[c(1:last_substance_index),])
    if (current_row > nrow(gc_peak_list[[s]])) {
        # stop when the last row was handled
        current_row <- 'stop'
    }
    }
    if (exists("pb")) close(pb)

    return(gc_peak_list)
}

