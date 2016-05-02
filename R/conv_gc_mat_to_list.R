#' transform GC data to list of individual GC matrices 
#' 
#' Transforms a Matrix containing GC data (retention times, peak area, peak hight etc) in adjacent
#' columns.
#' 
#' @param all_gc_mat matrix or data.frame containing GC data (retention times, peak area, peak hight etc) for all individuals
#'  in adjacent columns. The first column for all individuals has to be the retention time
#' @param ind_names character vector with names of individuals
#' @param var_names character vector with names of columns from all_gc_mat, 
#' such as c("RT", "Area", "Height"). The first name has to specify the retention time
#' 
#' @return 
#' List of data.frames. Each list element is the GC data for one individual 
#'
#' @references 
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      
#' @export
#' 


conv_gc_mat_to_list <- function(all_gc_mat, ind_names, var_names) {
    extract <- seq(from = 1, to = ncol(all_gc_mat), by = length(var_names))
    chromatograms <- lapply(extract, function(x) all_gc_mat[, x:(x+length(var_names)-1)])
    names(chromatograms) <- ind_names
    chromatograms <- lapply(chromatograms, rename_cols, var_names)
}