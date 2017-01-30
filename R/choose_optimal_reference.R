#' Picks the most suitable Reference for systematic shift corrections
#'
#' @description
#' Selects the chromatogram that shows the highest average similarity to all other samples a the reference for linear corrections of systematic shifts in peak retention times.
#'
#' @details
#' In order to correct systematic errors in peak retention times the most suitable reference chromatogram is selected based on the objective criterion of the highest average similarity to all other chromatograms in the data. Precisely, the median number of shared peaks between a sample and all other samples is estimated for any chromatogram in \code{gc_peak_list}. The most suitable referene is then defined as the sample showing the highest median similarity score.
#'
#' @param gc_peak_list
#' List of individual samples, where samples a data frames of numerical variables in columns. A column needs to contain retention times of peaks.
#'
#' @param rt_col_name
#' Name of the column in data frames of \code{gc_peak_list} that contains retention times.
#'
#' @return
#' returns a list with following elements
#' \item{sample}{Name of the sample with the highest average similarity to all other samples}
#' \item{score}{Median number of shared peaks with other samples}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import stats
#' @export
#'
choose_optimal_reference <- function(gc_peak_list = NULL, rt_col_name = NULL){
    if (is.null(rt_col_name)) stop("Specify retention time column")
    if (is.null(gc_peak_list)) stop("Provide a list of sampples")
    ## get the median scores for shared peaks
    x <- df_median_sim_score(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name)
    index <- which(x[["score"]] == max(x[["score"]]))
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

### define internal functions
## Scores for all samples
df_median_sim_score <- function(gc_peak_list,rt_col_name){
x <- data.frame(score = rep(NA,length(gc_peak_list)),sample = names(gc_peak_list))
    for(i in 1:length(gc_peak_list)) {
    x[["score"]][i] <- median_sim_score(gc_peak_list = gc_peak_list, ref_df = gc_peak_list[[i]],rt_col_name = rt_col_name)
    }
return(x)
    }

## Median similiarity to all other chromas
median_sim_score <- function(gc_peak_list,ref_df,rt_col_name){
    return(stats::median(unlist(lapply(gc_peak_list, chrom_sim_score,ref_df = ref_df,rt_col_name = rt_col_name))))
}

## Function comparin two samples
chrom_sim_score <- function(gc_peak_df, ref_df, rt_col_name,error=0) {
    ## get the rts of the reference
    ref_chroma <- ref_df[[rt_col_name]][!is.na(ref_df[[rt_col_name]])]
    sample_chroma <- gc_peak_df[[rt_col_name]][!is.na(gc_peak_df[[rt_col_name]])]
    ## Initialise vector tracking shared peaks
    peaks <- 0
    ## loop through the sample
    for (k in 1:length(sample_chroma)) {
        ## loop through the Reference Chromatogram
        for (l in 1:length(ref_chroma)) {
            temp_peak <- sample_chroma[k]
            ref_peak <- ref_chroma[l]
            ## Avoid comparison with cases of RT=0
            if (temp_peak != 0) {
                if ((temp_peak <= ref_peak + error) & (temp_peak >= ref_peak - error)) {
                    peaks <- peaks + 1
                }
            }
        }
    }
    return(peaks)
}
