#' Select the optimal reference for full alignments of peak lists
#'
#' @description
#' Full alignments of peak lists require the specification of a fixed reference to which all other samples are aligned to. This function provides an simple algorithm to find the most suitable sample among a dataset. The so defined reference can be used for full alignments using \code{\link{linear_transformation}}. The functions is evoked internally by \code{\link{align_chromatograms}} if no reference was specified by the user.
#'
#' @details
#' Every sample is considered in determining the optimal reference in comparison to all other samples by estimating the similarity to all other samples. For a reference-sample pair, the deviation in retention times between all reference peaks and the always nearest peak in the sample is summed up and divided by the number of reference peaks. The calculated value is a similarity score that converges to zero the more similar reference and sample are. For every potential reference, the median score of all pair-wise comparisons is used as a similarity proxy. The optimal sample is then defined by the minimum value among these scores. This functions is used internally in \code{\link{align_chromatograms}} to select a reference if non was specified by the user.
#'
#' @inheritParams align_chromatograms
#'
#' @examples
#' ## 1.) input is a list
#' ## using a list of samples
#' data("peak_data")
#' ## subset for faster processing
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
choose_optimal_reference <- function(data = NULL, rt_col_name = NULL, sep = "\t") {
    if (is.null(rt_col_name)) stop("Specify retention time column")
    if (is.null(data)) stop("Provide a list of sampples or the path to a text file")
    gc_peak_list <- data
    if (is.character(data)) gc_peak_list <- read_peak_list(data = data, rt_col_name = rt_col_name)

    # Currently only one method is supported
    method <- "Deviance"


    ## get the median scores for shared peaks
    x <- df_median_sim_score(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name, method = method)

    ## take the best, depending on the method chosen
    if (method == "Match") {
        index <- which(x[["score"]] == max(x[["score"]]))
    } else if (method == "Deviance") {
        #index <- which(min(x[["score"]]/x[["n_peaks"]]) == min(x[["score"]]/x[["n_peaks"]]))
        # Sun Jun 19 22:40:46 2022 ------------------------------
        # Bugfix thanks to hebertodelrio on GitHub!
        index <- which(x[["score"]]/x[["n_peaks"]] == min(x[["score"]]/x[["n_peaks"]]))
    }

    ## If more than one would get the same score, take the most central run
    if (length(index) > 1) {
    ## Odd or even number of samples determines the most central element
    centre <- ifelse(length(gc_peak_list) %% 2,length(gc_peak_list)/2 + 0.5,length(gc_peak_list)/2)
    diffs <- abs(centre - index)
    ## If still more than one index is equally well suited, i.e. Number 2 & 4
    index <- which(diffs == min(diffs))[1]
    }
    return(list(sample = as.character(x[["sample"]][[index]]), score = x[["score"]][[index]]/x[["n_peaks"]][[index]]))

    }
#### internal functions ####
#### ################## ####
df_median_sim_score <- function(gc_peak_list, rt_col_name, method) {

    pbapply::pboptions(char = "+", style = 1) # set up timer
    temp <- pbapply::pblapply(gc_peak_list, function(x) median_sim_score(gc_peak_list = gc_peak_list, ref_df = x, rt_col_name = rt_col_name, method = method))
    temp <- do.call("rbind", temp)
    ## number of peaks per sample
    temp_gc_peak_list <- remove_gaps(gc_peak_list = gc_peak_list, rt_col_name = rt_col_name)
    temp <- data.frame(score = as.vector(temp),
                       sample = rownames(temp),
                       n_peaks = as.vector(unlist(lapply(temp_gc_peak_list, nrow))))
    return(temp)
}

## Median similiarity to all other chromas
median_sim_score <- function(gc_peak_list, ref_df, rt_col_name, method = method){
    return(stats::median(unlist(lapply(gc_peak_list, chrom_sim_score,ref_df = ref_df,rt_col_name = rt_col_name, method = method))))
}

## Function comparing two samples
chrom_sim_score <- function(gc_peak_df, ref_df, rt_col_name, method = method) {
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
