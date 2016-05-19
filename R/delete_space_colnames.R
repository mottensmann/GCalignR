#' Remove spaces in column names
#' @description
#' If present, spaces within column names are eliminated
#'
#' @param data data.frame containing GC data.
#' @return
#' data.frame containing GC data with removed spaces in column names
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
delete_space_colnames <- function(data) {
    names(data) <- stringr::str_replace_all(names(data), " ", "")
    data
}

