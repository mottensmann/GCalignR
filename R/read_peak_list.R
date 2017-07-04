#' Read content of a text file
#'
#' @param data a text file
#'
#' @inheritParams align_chromatograms
#'
#' @description reads the content of text file that is formatted as described in \code{\link{align_chromatograms}} and creates a peak list.
#'
#' @return a list of peak lists for each sample
#'
#' @keywords internal
#'
#' @export
#'
read_peak_list <- function(data, sep, rt_col_name) {
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
    gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F)
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
