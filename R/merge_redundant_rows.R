#' Merge redundant rows
#'
#' @description
#' Sometimes, redundant rows (i.e. groups of resembling a homologous peak) remain in an aligned dataset. This is the case when two or more adjacent rows exhibit a difference in the mean retention time that is greater than \code{min_diff_peak2peak}, the parameter that determines a threshold below that redundancy is checked within \code{\link{align_chromatograms}}. Therefore, this function allows to raise the threshold for a post processing step that groups the homologous peaks together without the need of repeating a potentially time-consuming alignment with adjusted parameters.
#'
#' @details
#' Based on the value of parameter \code{threshold}, possibly redundant rows are identified by comparing mean retention times. Next, rows are checked for redundancy. When one or more samples contain peaks in a pair of compared rows, no redundancy is existent and the pair is skipped.
#'
#' @param data
#' An object of class "GCalign". See \code{\link{align_chromatograms}} for details.
#'
#' @param min_diff_peak2peak
#' A numerical giving a threshold in minutes below which rows of similar retention time are checked for redundancy.
#'
#' @return a list of two items
#' \item{GCalign}{input data with updated input to \code{\link{gc_heatmap}}}
#' \item{peak_list}{a list of data frames containing the updated dataset}
#'
#' @examples
#' ## Load example dataset
#' data("peak_data")
#' ## Subset for faster processing
#' peak_data <- peak_data[1:3]
#' peak_data <- lapply(peak_data, function(x) x[1:50,])
#' ## align data whith strict parameters
#' out <- align_chromatograms(peak_data, rt_col_name = "time",
#' max_diff_peak2mean = 0.01, min_diff_peak2peak = 0.02)
#' ## relax threshold to merge redundant rows
#' out2 <- merge_redundant_rows(data = out, min_diff_peak2peak = 0.05)
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
merge_redundant_rows <- function(data, min_diff_peak2peak = NULL) {
    if (!methods::is(data, "GCalign")) stop("Only data of type GCalign is supported")
    if (is.null(min_diff_peak2peak)) stop("Specify an numeric threshold value in minutes")
    gc_peak_list_aligned <- data[["aligned_list"]]

    rt_col_name <- data[["Logfile"]][["Call"]][["rt_col_name"]]

    ## add linear shifts again
    gc_peak_list_aligned <- add_linshifts2(dx = gc_peak_list_aligned, rt_col_name = rt_col_name, Logbook = data$Logfile)


    # number of peaks prior to executing this function
    N <- nrow(gc_peak_list_aligned[[1]])


    cat("Evaluate redundancy ...\n")
    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned,
                                                  min_diff_peak2peak = min_diff_peak2peak,
                                                  rt_col_name = rt_col_name)
    n <- N - nrow(gc_peak_list_aligned[[1]])
    cat(paste0(as.character(n)," rows were identified as redundant and were merged\n"))

    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)
    # update heatmap and aligned_list
    data[["heatmap_input"]][["aligned_rts"]] <-
        rt_extract(gc_peak_list = gc_peak_list_aligned,rt_col_name = rt_col_name)
    data[["aligned_list"]] <- gc_peak_list_aligned

    # update "aligned"
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))
    mean_per_row <- apply(rt_mat,1, function(x) if (all(x == 0)) 0 else mean(x[x != 0]))
    col_names <- names(gc_peak_list_aligned[[1]])
    output <- lapply(col_names, function(y) as.data.frame(do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[y]))))
    output <- lapply(output, function(x){
        names(x) <- names(gc_peak_list_aligned)
        x
    })

    output <- lapply(output, function(x){
        x <- cbind(mean_per_row, x)
        x
    })

    output <- lapply(output, function(x){
        names(x)[1] <- "mean_RT"
        x
    })
    names(output) <- col_names
    data[["aligned"]] <- output

    output <- list(GCalign = data, peak_list = gc_peak_list_aligned)

    return(output)
}







