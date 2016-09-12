#' align peaks individually among chromatograms
#'
#' @description
#' \code{align_peaks} allows to align similar peaks across samples so that
#' shared peaks are consistently located at the the same location (i.e. defined as the same substance).
#'  The order of chromatograms (i.e. data.frames in \code{gc_peak_list}) is randomized before each run
#'  of the alignment of algorithm. The main principle of this function is to reduce the variance in retention
#'  times within rows, thereby peaks of similar retention time are grouped together. Peaks that deviate
#'  signifantly from the mean retention times of the other samples are shifted to another row. At the begin
#'  of a row the first two samples are compared and seperated if required, then all other samples are included
#'  subsequently. If \code{n_iter>1} the whole alogorithm is repeated accordingly.
#'
#'@details
#'  For each row (i.e. peak) the algorithm compares the retention of every sample with those before. Starting
#'  with the second sample a comparison is done between the first and the second sample, then between the third
#'  and the two first ones etc. Whenever the current sample shows a deviation from the mean retention time of
#'  the previous samples a shift will either move this sample to the next row (i.e. retention time above average)
#'  or all other samples will be moved to the next row (i.e. retention time below average).
#'  If the retention time of the sample in evaluation shows no deviation within \code{-max_diff_peak2mean}:
#'  \code{max_diff_peak2mean} around the average retention time no shifting is done and the alogrithm takes the
#'  next sample.
#'
#' @param gc_peak_df
#' data.frame containing GC-data (e.g. retention time, peak area, peak height) of one sample. Variables are stored in columns.
#' Rows represent distinct peaks. Retention time is a required variable.
#'
#' @param gc_peak_list
#' List of data.frames. Each data.frame contains GC-data (e.g. retention time, peak area, peak height) of one sample. Variables are stored in columns.
#' Rows represent distinct peaks. Retention time is a required variable.
#'
#' @inheritParams align_chromatograms
#'
#' @param R
#' integer indicating the current iteration of the alignment step. Created by \link{align_chromatograms}.
#'
#' @return
#' a list of data.frames containing GC-data with aligned peaks.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @keywords internal
#'
#' @export
#'

align_peaks <- function(gc_peak_list, max_diff_peak2mean = 0.02, n_iter = 1, rt_col_name,R=1) {

    cat(paste('\n','\n','Iteration',as.character(R),'out of',as.character(n_iter),' ... '))
    shuffle_order <- sample(1:length(gc_peak_list)) # Generate random order
    gc_peak_list <- gc_peak_list[shuffle_order] # Shuffle

    current_row <- 1 # Start in Row 1
    while (current_row != 'Stop'){
# Loop through all rows, starting with one sample and adding others subsequentially
# Start with second sample
# Compare RT with mean of the others
# Shift current sample down, if RT is above mean
# Shift all previous samples, if RT is below
# Leave Samples with RT within error range at their positions

    for (S in 2:length(gc_peak_list)){
        current_RT <- gc_peak_list[[S]][current_row, rt_col_name] # Select a single retention time
        av_rt <- mean_retention_time_row(gc_peak_list, samples =1:(S-1), retention_row = current_row, rt_col_name)
    if(is.na(av_rt)){ # Selected peak is the only one in that row
        av_rt <- current_RT  # Do not shift
    }
    if(current_RT==0){ # Do not shift, if current has no substance in that row
        current_RT <- av_rt # Indicates not to shift anything
    }
    if (current_RT > av_rt + max_diff_peak2mean) { # Peak RT > Mean of Population
        gc_peak_list <- shift_rows(gc_peak_list,S,current_row)
    } else if (current_RT < (av_rt - max_diff_peak2mean)) { # Peak RT < Mean of Population
        for(J in 1:(S-1)){
            if(gc_peak_list[[J]][current_row, rt_col_name] <= (current_RT + max_diff_peak2mean)){
              # Do nothing, substance? position is valid
            } else {
                gc_peak_list <- shift_rows(gc_peak_list,J,current_row)
            }
        }
    }

}
    current_row <- current_row+1 # Next Peak to take
    gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)# Make all equal in length
    last_substance_index <- max(which(mean_retention_times(gc_peak_list, rt_col_name)==max(mean_retention_times(gc_peak_list, rt_col_name),na.rm=T)))

    gc_peak_list <- lapply(gc_peak_list, function(x) x[c(1:last_substance_index), ]) # Remove appended zeros
    if (current_row>nrow(gc_peak_list[[S]])){
        current_row <- 'Stop' # Signal the end of the iteration
    }
}

    return(gc_peak_list)
}


# Specifying internal functions
############################################################################################

similar_peaks <- function(average_rts, min_diff_peak2peak=0.05){

    difference <- rep(NA, (length(average_rts)-1)) # Estimate Difference between adjacent rows

    for (i in 2:length(average_rts)){
        difference[i] <- average_rts[i]-average_rts[i-1]
    }
    similar <- which(difference <= min_diff_peak2peak) # Which rows differ less the MinDistance ?
    return(similar)
}

mean_retention_times = function(gc_peak_list, rt_col_name){
    n_substance <- nrow(gc_peak_list[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance,
                         function(x) mean_retention_time_row(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    out

}

is_redundant <- function(redundant, criterion="strict"){
    # Indicates by a binary output variable (1/0) if rows should be merged
    # Methods: Strict: A single sample with two peaks prevents merging
    #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
    ToMerge <- 0
    if (criterion == "strict"){
        if(sum(redundant)/length(redundant) == 1){
            ToMerge <- 1
        }
    } else if (criterion == "proportional"){
        if(sum(redundant)/length(redundant)>=0.95){
            ToMerge <- 1
        }
    }
    ToMerge
}

check_redundancy <- function(gc_peak_df, similar_peaks, rt_col_name){
    # If only one of two neighbouring rows contain a substance
    # they are redundant, coded by a One
    row1 <- gc_peak_df[similar_peaks-1, rt_col_name] # Extract previous row
    row2 <- gc_peak_df[similar_peaks, rt_col_name] # Extract current row
    redundant <- 0
    if (row1==0 | row2==0){
        redundant <- 1
    }
    return(redundant)
}
merge_rows <- function(gc_peak_df, to_merge, criterion="strict", rt_col_name,conc_col_name){
    # Check always the row containing just zeros, in case of zeros in both, just delete one of them
    # To Merge == Last row of a similar pair
    Row1 <- to_merge-1 # Previous Row
    Row2 <- to_merge # Current Row
    R1 <- gc_peak_df[Row1, rt_col_name]
    R2 <- gc_peak_df[Row2, rt_col_name]
    if (criterion=="strict"){

        if (R1 == 0){
            #  Delete Row1, if no peak exists
            gc_peak_df <- gc_peak_df[-Row1,]
        } else if (R2 == 0){
            # Delete Row2
            gc_peak_df <- gc_peak_df[-Row2,]
        }
    }

    if(criterion=="proportional"){ # Take the larger of two peaks, defined by the are of the peaks

        if (gc_peak_df[Row1,conc_col_name] >= gc_peak_df[Row2,conc_col_name]){ # Keep Peak1, skip Peak2
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


shift_rows = function(chromatograms, current_sample_index, retention_row){
    # Shift Rows to sort substances by rows
    # List of Chromatograms
    # Index is the position of the Object in Chromatogram, which needs to be shifted
    n_col<- ncol(chromatograms[[1]])
    zeros <- as.data.frame(matrix(0,nrow=1,ncol=n_col))           # Create zeros
    colnames(zeros) <- names(chromatograms[[1]])

    chroma_temp <-  chromatograms[[current_sample_index]]

    if(retention_row != 1){
        chroma_temp<- rbind(chroma_temp[1:(retention_row-1), ], zeros,
                            chroma_temp[retention_row:nrow(chroma_temp), ])
    } else {
        chroma_temp <- rbind(zeros, chroma_temp)
    }
    # submit the new list, where the shift was applied to the Object at position current_sample_index
    chromatograms[[current_sample_index]] <- chroma_temp
    chromatograms
}


merge_redundant_peaks <- function(gc_peak_list,min_diff_peak2peak=0.05, rt_col_name,conc_col_name,criterion="strict"){

    merging <- 'Start'
    while(merging != 'Stop'){

        average_rts <- mean_retention_times(gc_peak_list, rt_col_name) # Average RTs, after merging rows
        similar <- similar_peaks(average_rts, min_diff_peak2peak)    # remaining similarities
        counter <- 1
        while (counter!='Stop'){ # allows to update after each merge
            if (length(similar) == 0){ # break loop if no similar peaks are remaining
                merging <- "Stop"
                break
            }
            redundant <- sapply(lapply(gc_peak_list, check_redundancy,similar_peaks=similar[counter], rt_col_name=rt_col_name), as.vector) #Check first position

            # checks whether all individuals just have substances in one of the rows ("strict")
            # checks whether at least 95% of individuals just have substances in one of the rows ("proportional")
            to_merge <- is_redundant(redundant = redundant, criterion=criterion)

            if (to_merge == 1){ # only merge if criterion proves redundancy of one of the rows
                gc_peak_list <- lapply(gc_peak_list, merge_rows, to_merge = similar[counter], criterion, rt_col_name,conc_col_name)
                # check2 <- similar[counter]
                counter <- 'Stop'
                # check <- do.call(cbind, lapply(gc_peak_list, function(x) x$RT))

            } else if  (to_merge == 0) {
                counter <- counter+1

                if (counter>length(similar)){
                    merging <- 'Stop'
                    counter <- 'Stop'
                }
            }
        }
    }

    return(gc_peak_list)
}

# End defining internal functions
##############################################################################################
