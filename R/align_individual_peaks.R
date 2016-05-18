#' algorithm to adjust unreliability of chromatograms on an individual peak basis
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param error_span maximum error under which two retention times are still counted as the same peak
#' @param n_iter number of iterations for the algorithm. One iteration is usually optimal, but
#'               for large datasets more iteration might be better. Has to be evaluated properly.
#' @param rt_col_name column name of retention time
#'
#' @return
#' \item{aligned chromatograms}{List of aligned peaks}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'

align_individual_peaks <- function(chromatograms, error_span = 0.02, n_iter = 1, rt_col_name) {

    for (R in 1:n_iter){
        cat(paste('Iteration',as.character(R),'out of',as.character(n_iter),'\n')) # Need to test whether it workds

        shuffle_order <- sample(1:length(chromatograms))
        chromatograms <- chromatograms[shuffle_order] # Shuffle
        # Give initial values to some variables
        # error_span <- 0.021 # allowed deviation around the mean of each row
        current_row <- 1

        while (current_row != 'Stop'){
            # Loop through all rows, starting with one sample and adding others subsequentially
            # Start with second sample
            # Compare RT with mean of the others
            # Shift current sample down, if RT is above mean
            # Shift all previous samples, if RT is below
            # Leave Samples with RT within error range at their positions

            for (S in 2:length(chromatograms)){
                current_RT <- chromatograms[[S]][current_row, rt_col_name]
                av_rt <- mean_of_samples(chromatograms, samples =1:(S-1), retention_row = current_row, rt_col_name)

                # if all rows are 0 ?
                if(is.na(av_rt)){
                    # Do not shift
                    av_rt <- current_RT
                }

                if(current_RT==0){
                    # Do not shift, if current has no substance in that row
                    current_RT <- av_rt
                }

                if (current_RT > av_rt + error_span) {
                    chromatograms <- shift_rows(chromatograms,S,current_row)

                } else if (current_RT < (av_rt - error_span)) {
                    for(J in 1:(S-1)){
                        if(chromatograms[[J]][current_row, rt_col_name] <= (current_RT + error_span)){
                            # Do nothing, substance?position is valid
                        } else {
                            chromatograms <- shift_rows(chromatograms,J,current_row)
                        }
                    }
                }

            }
            current_row <- current_row+1
            chromatograms <- lapply(chromatograms, matrix_append, chromatograms)# Make all equal in length
            last_substance_index <- max(which(mean_per_row(chromatograms, rt_col_name)==max(mean_per_row(chromatograms, rt_col_name),na.rm=T)))

            # rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x[, rt_col_name]))
            # max_rt_row <- max(rt_mat)
            #
            # last_substance_index <- lapply(chromatograms, function(x))
            # Remove tail of all-zero rows
            chromatograms <- lapply(chromatograms, function(x) x[c(1:last_substance_index), ]) # Remove appended zeros

            if (current_row>dim(chromatograms[[S]])[1]){
                current_row <- 'Stop' # Signal the end of the iteration
            }

        }

        # for evaluation of number of iteration, maybe create new function for this.
        # Length[R+1] <- max(unlist(lapply(chromatograms,function(x) out <- nrow(x))))
        # Variation[R+1] <- mean(var_per_row(chromatograms),na.rm = T)
    }

chromatograms
}

