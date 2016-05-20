#' Delete empty rows
#'
#' @description
#' \code{del_empty_rows} deletes all rows of data.frame that are empty. This is true for rows
#' with an average retention time of value NA (i.e. no peak is present in any data.frame) at this
#' location.
#'
#' @inheritParams matrix_append
#' @inheritParams similar_peaks
#'
#' @param average_rts numeric
#'@export

delete_empty_rows <- function(gc_peak_df, average_rts){
gc_peak_df <- gc_peak_df[!is.na(average_rts), ]
gc_peak_df
}
