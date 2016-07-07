#' cuts peaks below and above specified retention times
#'
#' @description By specifying cut-off thresholds retention times above (\code{highw}) or below (\code{high}) these peaks can
#'              be removed from the data in order to retain a subset of peaks belonging to
#'              specific retention times. All cut-off values are given in minutes.
#'
#'
#' @param rt_col_name
#' character string with the name of the retention time column in \code{gc_peak_df}
#'
#' @param low
#' lower threshold for retention times. RTs higher than \code{low} will be kept
#'
#' @param high
#' upper threshold for retention times. RTs lower than \code{high} will be kept
#'
#' @inheritParams matrix_append
#'
#' @return
#' data.frame without retention times inside \code{low}:\code{high}
#'
#' @details In addition to cutting retention times, this function removes rows of \code{gc_peak_df} containing
#'          purely NA. These are artefacts added by \code{\link[utils]{read.table}} when importing the data.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @keywords internal
#' @export




rt_cutoff <- function(gc_peak_df, rt_col_name, low=NULL, high=NULL){
  # RetentionCutoff removes all Retention Times below the Threshold specified by Low (default 8s).
  # In addition Retention Times above a time defined by the Value of High (Default is Null)
  # can be applied.
    highrow <- nrow(gc_peak_df)
    lowrow <- 1
    if (!is.null(low)){
        lowrow <- min(which(gc_peak_df[[rt_col_name]] > low))
    }
    if (!is.null(high)){
        highrow <- max(which(gc_peak_df[[rt_col_name]] < high))
    }

    gc_peak_df <- gc_peak_df[lowrow:highrow, ]
    out <- gc_peak_df[!is.na(gc_peak_df[, rt_col_name]), ]
    return(out)
}
