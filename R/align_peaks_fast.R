#'align peaks individually among chromatograms
#'
#'@description \strong{align_peaks_fast} is a fast implementation to sort peaks strictly based on the exact retention time values across samples. See \link{align_peaks} for details on the standard approach.
#'
#'@details Creates a matrix with the number of rows corresponding to the unique number of retention times across the dataset. Then, for each sample peaks are mapped to the corresponding row (=retention time) of the matrix and in parallel \code{gc_peak_list} is updated.
#'
#'@param gc_peak_list List of data.frames. Each data.frame contains GC-data
#'  (e.g. retention time, peak area, peak height) of one sample. Variables are
#'  stored in columns. Rows represent distinct peaks. Retention time is a
#'  required variable.
#'
#'@inheritParams align_chromatograms
#'
#'@return a list of data.frames containing GC-data with aligned peaks.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@keywords internal
#'
align_peaks_fast <- function(gc_peak_list, rt_col_name) {

    ## 1.) Get unique retention times
    ## -------------------------------------------------------------------------
    rt <- lapply(gc_peak_list, function(x) as.numeric(x[[rt_col_name]]))
    rt <- sort(unique(unlist(rt)))

    ## 2.) Create template for the output
    ## -------------------------------------------------------------------------
    df_out <- matrix(0, nrow = length(rt), ncol = ncol(gc_peak_list[[1]]))

    ## 2.) Find matches between mat and individual df in gc_peak_list
    ## -------------------------------------------------------------------------
    gc_peak_list <- lapply(gc_peak_list, function(df) {
        x <- match(df[[rt_col_name]], rt)
        df_out[x,] <- as.matrix(df)

        df_out <- as.data.frame(df_out)
        names(df_out) <- names(df)
        return(df_out)
    })
    return(gc_peak_list)
}


