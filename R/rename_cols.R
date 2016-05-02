#' Renames table columns in GC list
#' 
#' @param data \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time
#' @param var_names character vector with names of columns for \code{data}, such as c("RT", "Area",
#'   "Height"). The first name has to specify the retention time
#'   
#' @return \code{data.frame} with renamed columns
#' 
#' @references
#' 
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'   (meinolf.ottensmann@@web.de)
#'   
#' @export
#' 

rename_cols = function(data, var_names){
  names(data) <- var_names
  data
}


#' Deletes spaces when renaming table columns in GC list
#' 
#' @param data \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time
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
delete_space_colnames <- function(data) {
  names(data) <- stringr::str_replace_all(names(data), " ", "")
  data
}

