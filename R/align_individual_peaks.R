#' align peaks of chromatograms
#'
#' @description
#' \code{align_individual_peaks()} allows to align similar peaks across samples so that
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
#' @inheritParams matrix_append
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
#' @export
#'

align_individual_peaks <- function(gc_peak_list, max_diff_peak2mean = 0.02, n_iter = 1, rt_col_name,R=1) {

    #for (R in 1:n_iter){

        cat(paste('\n','\n','Iteration',as.character(R),'out of',as.character(n_iter),' ... ')) # Need to test whether it works
    set.seed(999) # Remove this after testing the algorithm !!!!!!!!!!!!!!
    shuffle_order <- sample(1:length(gc_peak_list))
        gc_peak_list <- gc_peak_list[shuffle_order] # Shuffle

        current_row <- 1

        while (current_row != 'Stop'){
            # Loop through all rows, starting with one sample and adding others subsequentially
            # Start with second sample
            # Compare RT with mean of the others
            # Shift current sample down, if RT is above mean
            # Shift all previous samples, if RT is below
            # Leave Samples with RT within error range at their positions

            for (S in 2:length(gc_peak_list)){
                current_RT <- gc_peak_list[[S]][current_row, rt_col_name]
                av_rt <- mean_retention_time_row(gc_peak_list, samples =1:(S-1), retention_row = current_row, rt_col_name)

                # if all rows are 0 ?
                if(is.na(av_rt)){
                    # Do not shift
                    av_rt <- current_RT
                }

                if(current_RT==0){
                    # Do not shift, if current has no substance in that row
                    current_RT <- av_rt
                }

                if (current_RT > av_rt + max_diff_peak2mean) {
                    gc_peak_list <- shift_rows(gc_peak_list,S,current_row)

                } else if (current_RT < (av_rt - max_diff_peak2mean)) {
                    for(J in 1:(S-1)){
                        if(gc_peak_list[[J]][current_row, rt_col_name] <= (current_RT + max_diff_peak2mean)){
                            # Do nothing, substance?position is valid
                        } else {
                            gc_peak_list <- shift_rows(gc_peak_list,J,current_row)
                        }
                    }
                }

            }
            current_row <- current_row+1
            gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)# Make all equal in length
            last_substance_index <- max(which(mean_retention_times(gc_peak_list, rt_col_name)==max(mean_retention_times(gc_peak_list, rt_col_name),na.rm=T)))

            # rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[, rt_col_name]))
            # max_rt_row <- max(rt_mat)
            #
            # last_substance_index <- lapply(gc_peak_list, function(x))
            # Remove tail of all-zero rows
            gc_peak_list <- lapply(gc_peak_list, function(x) x[c(1:last_substance_index), ]) # Remove appended zeros

            if (current_row>dim(gc_peak_list[[S]])[1]){
                current_row <- 'Stop' # Signal the end of the iteration
            }

        }


    #}
return(gc_peak_list)
}

