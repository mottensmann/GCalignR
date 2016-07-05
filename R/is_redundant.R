#' Checks if similar rows (i.e. substances) are redundant
#'
#' @description
#' Indicates by a binary output variable (1/0) which rows will be merged in the next step
#'
#' @param
#' redundant \code{vector} of integers coding redundancy as 1=Yes & 0=No for every sample at a given
#'          pair of similar rows.
#'
#' @inheritParams merge_rows
#'
#' @details
#'  Methods: "strict": A single sample with two peaks prevents merging
#           "proportional": Merging is acceptabel if only 5 % of samples show two peaks
#' @return
#' \item{ToMerge}{\code{integer} indicating whether will be merge (1) or not (0)}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @keywords internal


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
