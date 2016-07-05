#' Transforms GC-data to list
#'
#' @description
#' Transforms a data.frame containing GC-data (e.g. retention times, peak area, peak height etc.) to a list.
#' Data belonging to individual samples are concatenated horizontally in \code{gc_data}. Every samples needs to
#' have the same number of columns in consistent order. This means for example columns 1:3 belong to \code{Ind A}
#' and columns 4:6 to \code{Ind B} and so forth.
#'
#' @param gc_data
#' data.frame containing GC-data for all individuals in adjacent columns. Samples need to have the same
#' number of columns (i.e variables) in the same the order.
#'
#' @param ind_names
#' character vector with names of individuals following the same order as in \code{gc_data}.
#'
#' @param var_names
#' character vector with names of columns from gc_data. For example:
#' \code{var_names =c("RT", "Area", "Height")}.
#'
#' @return
#' A list of data.frames containing GC-data of individual samples
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @export
#'

conv_gc_mat_to_list <- function(gc_data, ind_names, var_names) {
    extract <- seq(from = 1, to = ncol(gc_data), by = length(var_names))
    chromatograms <- lapply(extract, function(x) gc_data[, x:(x+length(var_names)-1)])
    names(chromatograms) <- ind_names

    rename_cols = function(data, var_names){
        names(data) <- var_names
        data
    }

    chromatograms <- lapply(chromatograms, rename_cols, var_names)
    return(chromatograms)
}
