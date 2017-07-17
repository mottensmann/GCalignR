#' Pick the Optimal Reference for Full Alignments of Peak Lists
#'
#' @description
#' Full alignments of peak lists require the specification of a fixed reference to which all other samples are aligned to. This function provides an simple algorithm to find the most suitable sample among a dataset. The so defined reference can be used for full alignments using \code{\link{linear_transformation}}.
#'
#' @details
#' Every sample is considered in determining the optimal reference in comparison to all other samples by estimating the similarity to all other samples in pair-wise comparisons. For each peak in the reference a search for a matching counterpart in the second sample is conducted. Peaks are scored as shared among samples, when the retention times match within an user-defined error margin. Thereby, for each sample the median number of shared peaks with all other samples is estimated as similarity proxy. The optimal sample is then defined by the maximum value among these scores. This functions is used internally in \code{\link{align_chromatograms}} to select a reference if non was specified by the user.
#'
#' @param gc_peak_list
#' A list of data frames. Each data frame contains the peak list of a single sample with one column defining retention times (see \code{rt_col_name}).
#'
#' @param rt_col_name
#' A character giving the name of the column in data frames of that contains retention times.
#'
#' @param error
#' Allowed error in matching retention times. By default with \code{error = 0} retention times need to match with a precision of two decimals.
#'
#' @param method
#' Method used to estimate pair-wise similarity
#'
#' @examples
#' data("peak_data")
#' ## Subset for faster processing
#' peak_data <- peak_data[1:3]
#' choose_optimal_reference(peak_data, rt_col_name = "time")
#'
#' @return
#' A list with following elements
#' \item{sample}{Name of the sample with the highest average similarity to all other samples}
#' \item{score}{Median number of shared peaks with other samples}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import stats
#'
#' @export
#'
choose_optimal_reference <- function(gc_peak_list = NULL, rt_col_name = NULL, error = 0, method = c("Match","Deviance")){
    if (is.null(rt_col_name)) stop("Specify retention time column")
    if (is.null(gc_peak_list)) stop("Provide a list of sampples")
    method <- match.arg(method)


    ## get the median scores for shared peaks
    x <- df_median_sim_score(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name, error = error, method = method)

    ## take the best, depending on the method choose
    if (method == "Match") {
        index <- which(x[["score"]] == max(x[["score"]]))
    } else if (method == "Deviance") {
        index <- which(x[["score"]] == min(x[["score"]]))
    }

    ## If more than one would get the same score, take the most central run
    if (length(index) > 1) {
    ## Odd or even number of samples determines the most central element
    centre <- ifelse(length(gc_peak_list) %% 2,length(gc_peak_list)/2 + 0.5,length(gc_peak_list)/2)
    diffs <- abs(centre - index)
    ## If still more than one index is equally well suited, i.e. Number 2 & 4
    index <- which(diffs == min(diffs))[1]
    }
    return(list(sample = as.character(x[["sample"]][[index]]), score = x[["score"]][[index]]))

    }
#### internal functions ####
#### ################## ####
df_median_sim_score <- function(gc_peak_list, rt_col_name, error, method) {

    pbapply::pboptions(type = "timer", char = "+", style = 1) # set up timer
    temp <- pbapply::pblapply(gc_peak_list, function(x) median_sim_score(gc_peak_list = gc_peak_list, ref_df = x, rt_col_name = rt_col_name, error = error, method = method))
    temp <- do.call("rbind", temp)
    temp <- data.frame(score = as.vector(temp), sample = rownames(temp))
    return(temp)
}

## Median similiarity to all other chromas
median_sim_score <- function(gc_peak_list, ref_df, rt_col_name, error, method = method){
    return(stats::median(unlist(lapply(gc_peak_list, chrom_sim_score,ref_df = ref_df,rt_col_name = rt_col_name, error = error, method = method))))
}

## Function comparing two samples
chrom_sim_score <- function(gc_peak_df, ref_df, rt_col_name, error = error, method = method) {
    ## get the rts of the reference
    ref_chroma <- ref_df[[rt_col_name]][!is.na(ref_df[[rt_col_name]])]
    sample_chroma <- gc_peak_df[[rt_col_name]][!is.na(gc_peak_df[[rt_col_name]])]

    if (method == "Match") {
        peaks  <- sum(unlist(lapply(X = ref_chroma, FUN = function(fx) {
            ifelse(test = min(round(abs(fx - sample_chroma),2)) == 0, yes = 1, no = 0)
        })))

    } else if (method == "Deviance") {
        peaks <- sum(unlist(lapply(X = ref_chroma, FUN = function(fx) {
            min(abs(fx - sample_chroma))
        })))

        }# end method
    return(peaks)
}

#### ################## ####

# check method
# if (method == "Match1") {
# ## Initialise vector tracking shared peaks
# peaks <- 0
# ## loop through the sample
# for (k in 1:length(sample_chroma)) {
#     ## loop through the Reference Chromatogram
#     for (l in 1:length(ref_chroma)) {
#         temp_peak <- sample_chroma[k]
#         ref_peak <- ref_chroma[l]
#         ## Avoid comparison with cases of RT=0
#         if (temp_peak != 0) {
#             if ((round(temp_peak,2) <= round(ref_peak,2) + error) & (round(temp_peak,2) >= round(ref_peak,2) - error)) {
#                 peaks <- peaks + 1
#             }
#         }
#     }
# }
# }
############################
