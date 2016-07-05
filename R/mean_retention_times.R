#' calculates mean retention time per row
#'
#' @description \code{mean_retention_times} obtains the mean retention time for every substance currently listed in
#' data.frames of \code{gc_peak_list}. Calculations are performed by \code{\link{mean_retention_time_row}}
#'
#' @inheritParams align_chromatograms
#' @inheritParams matrix_append
#'
#' @return
#' vector with mean retention times per peak in \code{gc_peak_list}.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @keywords internal
#' @export
#'

mean_retention_times = function(gc_peak_list, rt_col_name){
    n_substance <- nrow(gc_peak_list[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance,
                         function(x) mean_retention_time_row(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    out

}
