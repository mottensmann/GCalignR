#' Read content of a text file and convert it to a list
#'
#' @param data
#' A text file containing a peak list. See \code{\link{align_chromatograms}} for details.
#'
#' @param check
#' logical
#'
#' @inheritParams align_chromatograms
#'
#' @description
#' Reads the content of text file that is formatted as described in \code{\link{align_chromatograms}} and converts it to a list.
#'
#' @return A list of data frames containing peak data for every sample of \code{data}.
#'
#' @examples
#' path <- system.file("extdata", "simulated_peak_data.txt", package = "GCalignR")
#' x <- read_peak_list(data = path, rt_col_name = "rt")
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
read_peak_list <- function(data, sep = "\t", rt_col_name, check = T) {
    if (check == TRUE) {
        check <- check_input(data = data, sep = sep, rt_col_name = rt_col_name, message = FALSE)
    } else {
        check <- TRUE
    }
    if (check == FALSE) stop("data is malformed. See check_input and fix issues raised there.")
    ind_names <- readr::read_lines(data, n_max = 1)
    ind_names <- unlist(stringr::str_split(string = ind_names,pattern = sep))
    ind_names <- ind_names[ind_names != ""]
    # Get Variable Names
    col_names <- readr::read_lines(data, n_max = 1, skip = 1)
    col_names <- unlist(stringr::str_split(string = col_names,pattern = sep))
    col_names <- col_names[col_names != ""]
    col_names <- stringr::str_trim(col_names)
    ind_names <- stringr::str_trim(ind_names)
    # Get Data
    gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F, fill = T)
    # Remove just NA-rows
    gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
    # Remove empty rows
    gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
    # Transform to data frame
    gc_data <-  as.data.frame(gc_data)
    # gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))
    # convert to list
    gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)
    # convert retention times to numeric
    fx <- function(x,rt_col_name) {
        x[[rt_col_name]] <- as.numeric(x[[rt_col_name]])
        return(x)
    }
    gc_peak_list <- lapply(gc_peak_list,FUN = fx,rt_col_name = rt_col_name)
    # remove filled gaps in samples
    gc_peak_list <- lapply(gc_peak_list, FUN = function(x) {
        if (any(is.na(rowSums(x)))) {
          p <- as.vector(which(is.na(rowSums(x))))
          x <- x[-p,]
        }
        if (any(rowSums(x) == 0)) {
            p <- as.vector(which(rowSums(x) == 0))
            x <- x[-p]
        }
        return(x)
    })

}
