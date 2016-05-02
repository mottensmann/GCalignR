#' cuts retention time below and above specified times
#' 
#' @param data \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time
#' @param rt_col_name character string for name of retention time in gc table
#' @param low Lower threshold for retention times. RTs higher than \code{low} will be kept
#' @param high Upper threshold for retention times. RTs lower than \code{high} will be kept
#' 
#' @return 
#' gc table with retention times cut
#'
#' @references 
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      
#' @export
#' 




rt_cutoff <- function(data, rt_col_name, low=NULL, high=NULL){
  # RetentionCutoff removes all Retention Times below the Threshold specified by Low (default 8s).
  # In addition Retention Times above a time defined by the Value of High (Default is Null)
  # can be applied.
    highrow <- nrow(data)
    lowrow <- 1
    if (is.null(low) | low < 0){
        low <- 0
    } else {
        lowrow <- min(which(data[[rt_col_name]] > low))
    }
    if (is.null(high)| high < 0){
        high <- nrow(data)
    } else {
        highrow <- max(which(data[[rt_col_name]] < high))
    }
    
    out <- data[lowrow:highrow, ]

}
