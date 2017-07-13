#'align peaks individually among chromatograms
#'
#'@description \strong{align_peaks} allows to align similar peaks across samples
#'so that shared peaks are consistently located at the the same location (i.e.
#'defined as the same substance). The order of chromatograms (i.e. data.frames
#'in \code{gc_peak_list}) is randomized before each run of the alignment of
#'algorithm. The main principle of this function is to reduce the variance in
#'retention times within rows, thereby peaks of similar retention time are
#'grouped together. Peaks that deviate significantly from the mean retention times
#'of the other samples are shifted to another row. At the start of a row the
#'first two samples are compared and separated if required, then all other
#'samples are included consecutively. If \code{iterations > 1} the whole
#'alogorithm is repeated accordingly.
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
#'@param gc_peak_df data.frame containing GC-data (e.g. retention time, peak
#'  area, peak height) of one sample. Variables are stored in columns, rows
#'  represent peaks.
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
align_peaks <- function(gc_peak_list, max_diff_peak2mean = 0.02, iterations = 1, rt_col_name,R = 1) {
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
        pb <- txtProgressBar(min = 0, max = total, style = 3, char = "+", width = 80)
        #Sys.sleep(1)
        setTxtProgressBar(pb, current_row)
        ##############


        # Randomize the order of samples
        shuffle_order <- sample(1:length(gc_peak_list))
        gc_peak_list <- gc_peak_list[shuffle_order]

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
    close(pb) # end of the progress bars

    return(gc_peak_list)
}

# Specifying internal functions
############################################################################################

# Find similar rows
similar_peaks <- function(average_rts, min_diff_peak2peak = 0.05){

    difference <- rep(NA, (length(average_rts) - 1))
    for (i in 2:length(average_rts)) {
        difference[i] <- average_rts[i] - average_rts[i - 1]
    }
    # Find rows that show similar mean retention times
    similar <- which(difference <= min_diff_peak2peak)
    return(similar)
}

# calculate mean retention times
mean_retention_times <- function(gc_peak_list, rt_col_name) {
    n_substance <- nrow(gc_peak_list[[1]])
    out <- unlist(lapply(1:n_substance,
                         function(x) mean_retention_time_row(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    return(out)
}

# evaluate if a row of a pair of rows is redundant
is_redundant <- function(redundant, criterion="strict"){
    # Indicates by a binary output variable (1/0) if rows should be merged
    # Methods: Strict: A single sample with two peaks prevents merging
    #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
    ToMerge <- 0
    if (criterion == "strict") {
        if (sum(redundant)/length(redundant) == 1) {
            ToMerge <- 1
        }
    } else if (criterion == "proportional") {
        if (sum(redundant)/length(redundant) >= 0.95) {
            ToMerge <- 1
        }
    }
    ToMerge
}

# Just if rows can be merged
check_redundancy <- function(gc_peak_df, similar_peaks, rt_col_name){
    # If only one of two neighbouring rows contain a substance
    # they are redundant, coded with 1
    row1 <- gc_peak_df[similar_peaks - 1, rt_col_name]
    row2 <- gc_peak_df[similar_peaks, rt_col_name]
    redundant <- 0
    if (row1 == 0 | row2 == 0) {
        redundant <- 1
    }
    return(redundant)
}
# Apply the merge
merge_rows <- function(gc_peak_df, to_merge, criterion="strict", rt_col_name,conc_col_name){
    # Check always the row containing just zeros, in case of zeros in both, just delete one of them
    # To Merge == Last row of a similar pair
    Row1 <- to_merge - 1
    Row2 <- to_merge
    R1 <- gc_peak_df[Row1, rt_col_name]
    R2 <- gc_peak_df[Row2, rt_col_name]
    if (criterion == "strict") {
        if (R1 == 0) {
            #  Delete Row1, if no peak exists
            gc_peak_df <- gc_peak_df[-Row1,]
        } else if (R2 == 0) {
            # Delete Row2
            gc_peak_df <- gc_peak_df[-Row2,]
        }
    }

    if (criterion == "proportional") {
        # keep the peak with higher area
        if (gc_peak_df[Row1,conc_col_name] >= gc_peak_df[Row2,conc_col_name]) {
            gc_peak_df <- gc_peak_df[-Row2,]
        } else if (gc_peak_df[Row1,conc_col_name] < gc_peak_df[Row2,conc_col_name]) {
            gc_peak_df <- gc_peak_df[-Row1,]
        }
    }
    return(gc_peak_df)
}

mean_retention_time_row <- function(gc_peak_list, samples, retention_row, rt_col_name){
   xy <- function(x, retention_row, rt_col_name) {
       out <- x[retention_row, rt_col_name]
   return(out)
       }
    rts <- unlist(lapply(gc_peak_list[samples], xy,retention_row, rt_col_name))
    mean_rt <- mean(rts[!(rts == 0)], na.rm = TRUE)
return(mean_rt
       )
}

shift_rows = function(chromatograms, current_sample_index, retention_row) {
    n_col <- ncol(chromatograms[[1]])
    zeros <- as.data.frame(matrix(0,nrow = 1,ncol = n_col))
    colnames(zeros) <- names(chromatograms[[1]])
    chroma_temp <-  chromatograms[[current_sample_index]]

    if (retention_row != 1) {
        chroma_temp <- rbind(chroma_temp[1:(retention_row - 1),], zeros,
                            chroma_temp[retention_row:nrow(chroma_temp), ])
    } else {
        chroma_temp <- rbind(zeros, chroma_temp)
    }

    chromatograms[[current_sample_index]] <- chroma_temp
    return(chromatograms)
}


merge_redundant_peaks <- function(gc_peak_list,min_diff_peak2peak=0.05, rt_col_name, conc_col_name = NULL, criterion="strict"){

    merging <- 'start'
    while (merging != 'stop') {

        # calculate mean retention times
        average_rts <- mean_retention_times(gc_peak_list, rt_col_name)
        # update similarity assessment
        similar <- similar_peaks(average_rts, min_diff_peak2peak)
        counter <- 1

       while (counter != 'stop') {

           ##############
           ### timer ####
           ##############
           total <- length(similar)
           # create progress bar
           pb <- txtProgressBar(min = 0, max = total, style = 3, char = "+", width = 80)
           #Sys.sleep(1)
           setTxtProgressBar(pb, counter)
           ##############

            # stop when there are no redundancies
            if (length(similar) == 0) {
                merging <- "stop"
                break
            }
            redundant <- sapply(lapply(gc_peak_list, check_redundancy,similar_peaks = similar[counter], rt_col_name = rt_col_name), as.vector)

            to_merge <- is_redundant(redundant = redundant, criterion = criterion)
# prove if rows are mergeable
            if (to_merge == 1) {
                gc_peak_list <- lapply(gc_peak_list, merge_rows, to_merge = similar[counter], criterion, rt_col_name,conc_col_name)
                counter <- 'stop'
            } else if  (to_merge == 0) {
                counter <- counter + 1
                if (counter > length(similar)) {
                    merging <- 'stop'
                    counter <- 'stop'
                }
            }
       }
    }
    close(pb)
    return(gc_peak_list)
}
# End defining internal functions
##############################################################################################
