#' Remove spaces in column names
#' @description
#' If present, spaces within column names are eliminated
#'
#' @inheritParams conv_gc_mat_to_list
#'
#' @return
#' data.frame containing GC data with removed spaces in column names
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
delete_space_colnames <- function(gc_data) {
    names(gc_data) <- stringr::str_replace_all(names(gc_data), " ", "")
    gc_data
}

