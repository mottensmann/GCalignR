#' Check correct formatting of data input
#'
#'@description
#' Checks conformity between the input format and the requirements of GCalignR. Supported are
#' a \code{.txt} files or lists of data.frames. See \code{\link{align_chromatograms}} for details.
#'
#'@param data
#'       path to a data file or the name of a list in the Global Environment.
#'@param sep
#'The field separator character. Values on each line of the file are separated by this
#'character. The default is tab seperated (sep = '\\t'). See \code{sep} argument in \code{\link[utils]{read.table}} for details.
#'
#'@param ...
#'optional arguments used internally in GCalignR. See source code for details.
#'
#'@return TRUE if data is formatted correctly, warning and explanation if not.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@import magrittr stringr
#'
#' @examples
#' data(gc_peak_data)
#' gc_peak_data <- gc_peak_data[1:4]
#' check_input(gc_peak_data)
#'
#' @export
#'

check_input <- function(data, sep = "\t",...) {

    opt <- list(...) # optional parameters

    if (is.character(data)) { # Check if data is the path to a txt.file
        if (!stringr::str_detect(data, ".txt")) {
            stop("Data is not of the expected format. Specify a valid path to a .txt-file")
        }
        ind_names <- readr::read_lines(data, n_max = 1) %>% # Sample Names
            stringr::str_split(pattern = sep) %>%
            unlist()
        ind_names <- ind_names[ind_names != ""]

        col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>% # Variable Names
            stringr::str_split(pattern = sep) %>%
            unlist()
        col_names <- col_names[col_names != ""]
        col_names <- stringr::str_trim(col_names)
        ind_names <- stringr::str_trim(ind_names)

        gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F) # Peak Data
        gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ] # Remove just NA-rows
        gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))] # Remove empty rows
        gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))# Transform variables to numeric

        ################################
        # Check input for completeness #
        ################################
        if (!((ncol(gc_data) / length(col_names)) %% 1) == 0) stop("Number of data columns is not a multiple of the column names provided")
        if (!((ncol(gc_data) / length(col_names))  == length(ind_names))) stop("Number of sample names provided does not fit to the number of columns in the data")
        if(any(duplicated(ind_names)))warning("Avoid duplicates in sample names")
        gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names) # convert to list


    } else if (is.list(data)) { # If data is a list of data.frames
        if (!(all(unlist(lapply(data, is.data.frame))))) stop("Every Sample has to be a data.frame") # check that every element in the list is a data.frame
        if ((is.null(names(data)))) stop("Every data.frame needs to be named with the sample id")     # check all data.frames are named
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        if(any(table(col_names) != length(data))) stop("Every sample needs to have the same number of columns")
        col_names <- names(data[[1]])
        ind_names <- names(data)
        if(any(duplicated(ind_names)))warning("Avoid duplicates in sample names")
        gc_peak_list <- data
    }

    ##############################
    # Some checks further checks #
    ##############################
    if(any(names(opt)=="write_output")){
        if(any(!(opt[["write_output"]]%in%col_names))) stop("Strings in write_output need to indicate a variable contained in data!")
    }
    if(any(stringr::str_detect(string = ind_names, pattern = " "))) warning("Avoid whitespaces in Sample Names!")
    if(any(stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Sample Names should only contain Letters, Numbers and '_' ")
    if(any(stringr::str_detect(string = col_names, pattern = " "))) warning("Avoid whitespaces in Variable Names!")
    if(any(stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Variable Names should only contain Letters, Numbers and '_' ")

    format_error <- function(x){
        check_var_count <- function(x){
            mat <- as.matrix(x)
            L<-length(unique(colSums(!is.na(x)))==1)
            return(L)
        }
        y <- unlist(lapply(gc_peak_list,check_var_count))
        if(length(which(y!=1))>0){ #1
            out <- names(y[which(y!=1)])
            warning(paste(out,collapse = "; ") ," violate(s) the requirements.",call. = FALSE)
            stop("Every sample needs to have the same number of values for each variable!",call. = FALSE)
        }
    }

    format_error(gc_peak_list) # Checks that every sample has the same number of values per column
    cat("All checks passed!\nReady for processing with align_chromatograms")
    return("Passed")
}



