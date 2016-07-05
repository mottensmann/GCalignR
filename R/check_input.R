#' Check correct formatting of data input
#'
#'@description
#' Checks whether the data is formatted correctly. The data can either be a .txt file
#' or a list of data.frames. See \code{\link{align_chromatograms}}
#'
#'
#'@param data path to data file or list of data.frames
#'@param sep
#'The field separator character. Values on each line of the file are separated by this
#'character. The default is tab seperated (sep = '\\t'). See \code{sep} argument in \code{\link[utils]{read.table}} for details.
#'
#'
#'
#'@return TRUE if data is formatted correctly, warning and explanation if not.
#'
#'
#'  @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
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

check_input <- function(data, sep = "\t") {

    if (is.character(data)) {
        if (stringr::str_detect(data, ".txt")) {

            # extract names
            ind_names <- readr::read_lines(data, n_max = 1) %>%
                stringr::str_split(pattern = sep) %>%
                unlist()
            ind_names <- ind_names[ind_names != ""]    #.[. != ""]

            # some checks for the names
            if(any(stringr::str_detect(string = ind_names, pattern = " "))) warning("sample names should not contain whitespaces")
            if(any(stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("sample names should just contain letters, numbers and '_' ")
            if(any(duplicated(ind_names))) warning("sample names have to be unique")


            ######################
            # extract column names
            ######################
            col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>%
                stringr::str_split(pattern = sep) %>%
                unlist()

            col_names <- col_names[col_names != ""]    #.[. != ""]
            # some checks for the names
            if(any(stringr::str_detect(string = col_names, pattern = " "))) warning("column names should not contain whitespaces")
            if(any(stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("column names should just contain letters, numbers and '_' ")

            ########################################
            # remove leading and tailing whitespaces
            ########################################
            col_names <- stringr::str_trim(col_names)
            ind_names <- stringr::str_trim(ind_names)


            ##############
            # extract data
            ##############
            gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F)

            #####################
            # remove pure NA rows
            #####################
            gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]

            ############################################################
            # remove pure NA columns
            # Note: Putatively created when wrangling the data in excel
            # It is not visible, that empty columns are in the data!
            ############################################################
            gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
            ######################
            # transform to numeric
            ######################
            gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))

            #########################################
            # Check input for completeness and format
            #########################################
            # check 1
            if (!((ncol(gc_data) / length(col_names)) %% 1) == 0) warning("Number of data columns is not a multiple of the column names provided")
            # check 2
            if (!((ncol(gc_data) / length(col_names))  == length(ind_names))) warning("Number of sample names provided does not fit to the number of columns in the data")

            # check that each individual data has the same value count in each column
            check_val_counts <- function(x) {
                ifelse(!(length(unique(colSums(is.na(x)))) == 1), TRUE, FALSE)
            }

            format_error <- which(sapply(seq(from = 1, to = ncol(gc_data), by = length(col_names)),
                function(x) check_val_counts(gc_data[x:(x + length(col_names) - 1)])))
            if (length(format_error) != 0) warning(paste("samples", ind_names[format_error], "have an unequal number of values in their variables", sep = " "))
            warning("columns within samples have to have the same length")
        }

    } else if (is.list(data)) {
        # check if data is list of data.frames
        # check for data.frames
        if (!(all(unlist(lapply(data, is.data.frame))))) warning("Data object has to be a list, whereby each element is a data.frame with the GC peak data for an individual")
        # check whether all data.frames have a name
        if ((is.null(names(data)))) warning("Data object has to be a list, whereby each element has to be named with the ID of the respective individual")
        # check whether all data.frames have a unique name
        if(any(duplicated(names(data)))) warning("sample (data.frame) names have to be unique")
        # check whether all data.frames contain the same column names
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        if(any(table(col_names) != length(data))) warning("Each data.frame in the list has to have the same variable names (i.e. 'RT' 'area')")

        # check that each individual data has the same value count in each column
        check_val_counts <- function(x) {
            ifelse(!(length(unique(colSums(is.na(x)))) == 1), TRUE, FALSE)
        }

        format_error <- which(unlist(lapply(data, check_val_counts)))
        if (length(format_error) != 0) warning(paste("samples", ind_names[format_error], "have an unequal number of values in their variables", sep = " "))
    }

out <- TRUE
}



