#' cuts retention time below and above specified times
#'
#' @description By specifying cut-off thresholds retention times above (\code{highw}) or below (\code{high}) these values can
#'              be removed from the data in order to retain a subset of substances belonging to
#'              specific retention times. All cut-off values are given in minutes.
#'
#' @param data \code{data.frame} containing GC data (e.g. retention time, peak are, peak height etc.) for
#'   one individual. Variables are stored in adjacent columns.
#'
#' @param rt_col_name \code{character} string with the name of the retention time column
#'
#' @param low lower threshold for retention times. RTs higher than \code{low} will be kept
#'
#' @param high upper threshold for retention times. RTs lower than \code{high} will be kept
#'
#' @return
#' \item{data_cutted}{\code{data.frame} containing subset of the input data falling within the interval
#' defined by \code{low} and \code{high} respectively}
#'
#' @details In addition to cutting retention times, this function removes rows of \code{data} containing
#'          purely NA. These are artefacts added by \code{\link[utils]{read.table}}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export




rt_cutoff <- function(data, rt_col_name, low=NULL, high=NULL){
  # RetentionCutoff removes all Retention Times below the Threshold specified by Low (default 8s).
  # In addition Retention Times above a time defined by the Value of High (Default is Null)
  # can be applied.
    highrow <- nrow(data)
    lowrow <- 1
    if (!is.null(low)){
        lowrow <- min(which(data[[rt_col_name]] > low))
    }
    if (!is.null(high)){
        highrow <- max(which(data[[rt_col_name]] < high))
    }

    data <- data[lowrow:highrow, ]
    out <- data[!is.na(data[, rt_col_name]), ]
    return(out)
}
