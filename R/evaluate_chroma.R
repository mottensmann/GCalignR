#' calculates mean of retention times up to a certain sample (not including 0s)
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param samples indices of samples up to sample of interest (1:sample-1)
#' @param retention_row current retention time row to be compared
#' @param rt_col_name name of retention time column
#'
#' @return
#' mean of rts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'

mean_of_samples <- function(chromatograms, samples, retention_row, rt_col_name){
    rts <- unlist(lapply(chromatograms[samples], function(x) x[retention_row, rt_col_name]))
    mean_rt <- mean(rts[!(rts == 0)], na.rm = TRUE)
    ## round ?

}

#' calculates mean of retention times up to a certain sample (not including 0s)
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param samples indices of samples up to sample of interest (1:sample-1)
#' @param retention_row current retention time row to be compared
#' @param rt_col_name name of retention time column
#'
#' @return
#' var of rts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
var_of_samples <- function(chromatograms, samples, retention_row, rt_col_name){
    # Estimate the Variation in Retention Times within Rows
    # NA indicates that only one Substance exists
    rts <- unlist(lapply(chromatograms[samples], function(x) x[retention_row, rt_col_name]))
    mean_rt <- var(rts[!(rts == 0)], na.rm = TRUE)
}


#' calculates mean retention time per row
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @return
#' vector with mean retention times per row in chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export

mean_per_row = function(chromatograms, rt_col_name){
    n_substance <- nrow(chromatograms[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance,
                    function(x) mean_of_samples(chromatograms, 1:length(chromatograms), x, rt_col_name)))
    out

}



#' calculates mean retention time per row
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @return
#' vector with variance of retention times per row in chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)


var_per_row = function(chromatograms, rt_col_name){
    n_substance <- nrow(chromatograms[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance,
        function(x) var_of_samples(chromatograms, 1:length(chromatograms), x, rt_col_name)))
    out
}






