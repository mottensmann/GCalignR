#' Calculate mean retention for one peak
#'
#' @description
#' \code{mean_retention_time_row} calculates the average retention time of a row. Samples in \code{gc_peak_list}
#' lacking a peak (i.e retention time = 0) are not considered.
#'
#' @inheritParams matrix_append
#' @inheritParams align_chromatograms
#'
#' @param samples
#' indices of samples to consider in the calculation of the mean retention time.
#'
#' @param retention_row
#' index of the current retention time to be evaluate.
#'
#'
#' @return
#' mean retention time for peaks at index \code{retention_row}.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'

mean_retention_time_row <- function(gc_peak_list, samples, retention_row, rt_col_name){
    rts <- unlist(lapply(gc_peak_list[samples], function(x) x[retention_row, rt_col_name]))
    mean_rt <- mean(rts[!(rts == 0)], na.rm = TRUE)

}
