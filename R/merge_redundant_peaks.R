#' Merges redundant rows
#'
#'@description
#' \code{merge_redundant_peaks()} allows to merge neighbouring peaks (i.e substances) of similiar retention time. All peaks
#' that are similar in retention times (definded by \code{min_diff_peak2peak}) are combined, in case that
#' samples contain either one or none of two compared peaks, but never both. An optional procedure
#' allows to merge close peaks if only a small subset of samples contains both of them.
#'
#' @details
#' Peaks that are potentially redundant are selected by comparing the similarity of retention times with
#' respect to \code{min_diff_peak2peak}. Candidates are peaks that are closer to each other than this
#' threshold. After one pair of rows was merged, average retention times are calculated again. The algorithm
#' stops if no pair of similar peaks remains or none of those appears redundant (i.e. at least one sample
#' shows both peaks). As a last step, peaks that were not proved to be redundant by this strict criterion,
#' are merged if only a small proportion of samples (less than 5%) shows to conflicting peaks. In this
#' case the peak with the larger concentration (peak area or peak height, depends on \code{conc_col_name})
#' is retained and the other one deleted.
#'
#' @inheritParams matrix_append
#'
#' @inheritParams align_chromatograms
#'
#' @inheritParams merge_rows
#'
#' @return
#' a list of data.frames for individual samples with redundant peaks merged.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
merge_redundant_peaks <- function(gc_peak_list,min_diff_peak2peak=0.05, rt_col_name,conc_col_name,criterion="strict"){

    merging <- 'Start'
    while(merging != 'Stop'){

        average_rts <- mean_retention_times(gc_peak_list, rt_col_name) # Average RTs, after merging rows
        similar <- similar_peaks(average_rts, min_diff_peak2peak)    # remaining similarities

        counter <- 1

        # cat(length(similar),counter,"\n")
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
