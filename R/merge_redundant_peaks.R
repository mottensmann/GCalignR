#' Merges redundant rows
#'
#'@description
#' \code{merge_redundant_peaks()} allows to merge neighbouring peaks (i.e substances) of similiar retention time. All peaks
#' that are similar in retention times (definded by \code{min_diff_peak2peak}) are combined, in case that
#' samples contain either one or none of two compared peaks, but never both. An optional procedure
#' allows to merge close peak if only a small subset of samples contains both of them.
#'
#' @details
#' Peaks that are potentially redundant are determined by comparing the similarity of retention times with
#' respect to \code{min_diff_peak2peak}. Candidates are all peaks that are closer to each other than this
#' threshold. After one pair of rows was merged average retention times are calculated again. The algorithm
#' stops if no pair of similar peaks remains or none of those appears redunant (i.e. at least one sample
#' shows both peaks).
#'
#'
#' @inheritParams matrix_append
#'
#' @inheritParams align_chromatograms
#'
#' @return
#' a list of data.frames for individual samples with redundant peaks merged.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
merge_redundant_peaks <- function(gc_peak_list,min_diff_peak2peak=0.05, rt_col_name){
# removed average_rts from function input, since it is calculated within

    merging <- 'Start'
    while(merging != 'Stop'){


        # Updating merging criterions
        average_rts <- mean_retention_times(gc_peak_list, rt_col_name) # Average RTs, after merging rows
        similar <- similar_peaks(average_rts, min_diff_peak2peak)    # remaining similarities

        counter <- 1
        while (counter!='Stop'){
            # loop through vector of similar rows, until one shift was done

            # check if any chromatogram contains substances in both rows
            if (length(similar) == 0){
                merging <- "Stop"
                break
            }
            redundant <- sapply(lapply(gc_peak_list, check_redundancy, similar[counter], rt_col_name), as.vector) #Check first position

            # checks whether all individuals just have substances in one of the rows ("Strict")
            # checks whether at least 95% of individuals just have substances in one of the rows ("Proportional")
            criterion <- is_redundant(redundant = redundant, criterion="strict")

            if (criterion == 1){ # only merge if criterion proves redundancy of one of the rows
                gc_peak_list <- lapply(gc_peak_list, merge_rows, to_merge = similar[counter], criterion="zero", rt_col_name)

                # check2 <- similar[counter]
                counter <- 'Stop'

                # check <- do.call(cbind, lapply(gc_peak_list, function(x) x$RT))

            } else if  (criterion == 0) {
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
