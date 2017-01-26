#' Check input prior to processing in GCalignR
#'
#'@description
#' Checks conformity between the input format and the requirements of GCalignR. Data is accepted in form of the path to a text file (i.e. "data.txt") or a list of data frames. See \code{\link{align_chromatograms}} for details.
#'
#' @details
#' The data needs to fullfil certain formatting standards. Sample names are expected to contain just letters, numbers and underscores. Each sample hast to contain the same number of columns that is determined given the number of variables. It is mandatory that the values of all variables are numeric, i.e. they are allowed to contain numbers from 0-9 and "." as the only decimal character.
#'
#'@param data
#'       path to a data file or the name of a list in the global environment.
#'
#'@param plot
#'logical, if TRUE the distribution of peak numbers is plotted. Default is FALSE.
#'
#'@inheritParams align_chromatograms
#'
#'@param ...
#'optional arguments passed to methods, see \code{\link[graphics]{barplot}}. Reguires \code{plot == TRUE}.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@import magrittr stringr
#'
#'@return
#' If \code{plot = TRUE} a data frame containing sample names and the corresponding number of peaks is returned
#'
#' @examples
#' ## Checks format
#' check_input(peak_data)
#' ## Includes a barplot of peak numbers in the raw data
#' check_input(peak_data, plot = TRUE)
#'
#' @export
#'
check_input <- function(data,plot = FALSE, sep = "\t", ...) {

    # Preallocate a flag for passing the test. Every severe issues sets pass to FALSE
    pass = TRUE
    ## Get the name of "data" and optional parameters
    mcall = as.list(match.call())[-1L]
    opt <- list(...)

    ## Check files
    ## Check if data is the path to a txt.file
    if (is.character(data)) {
        if (!stringr::str_detect(data, ".txt")) {
            stop("Data is not of the expected format. Specify a valid path to a .txt-file")
        }
        ## Sample Names
        ind_names <- readr::read_lines(data, n_max = 1) %>%
            stringr::str_split(pattern = sep) %>%
            unlist()
        ind_names <- ind_names[ind_names != ""]
        ## Variable Names
        col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>%
            stringr::str_split(pattern = sep) %>%
            unlist()
        col_names <- col_names[col_names != ""]
        col_names <- stringr::str_trim(col_names)
        ind_names <- stringr::str_trim(ind_names)
        ## Get Peak Data
        gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F)
        ## Check that all values are numeric with proper decimal operator
        x <- as.vector(unlist(gc_data))
        na <- length(x[is.na(x)])
        y <- as.vector(as.data.frame(apply(gc_data, 2, as.numeric)))
        y <- as.vector(unlist(y))
        if (suppressWarnings(length(y[is.na(y)]) - na > 0)) stop("Data contains non-numeric data. Make sure not to use >,< as a decimal operator")
        ## Remove just NA-rows
        gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
        ## Remove empty rows
        gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
        ## Transform variables to numeric
        gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))

        ### Check input for completeness ###

        if (!((ncol(gc_data) / length(col_names)) %% 1) == 0) {
            pass <- FALSE
        warning("Number of data columns is not a multiple of the column names provided")
        }
        if (!((ncol(gc_data) / length(col_names))  == length(ind_names))) {
            pass <- FALSE
            warning("Number of sample names provided does not fit to the number of columns in the data")
        }
        if (any(duplicated(ind_names))) warning("Avoid duplicates in sample names")
        ## convert to list
        gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)

## If data is a list of data.frames
    } else if (is.list(data)) {
        ## check that every element in the list is a data.frame
        if (!(all(unlist(lapply(data, is.data.frame))))) {
            pass <- FALSE
            warning("Every Sample has to be a data.frame")
        }
        ## Adjust sizes of data frames, each has the same number of rows
        data <- lapply(data,matrix_append,gc_peak_list = data, val = "NA")
        ## Check that all values are numeric with proper decimal operator
        x <- as.vector(unlist(data))
        na <- length(x[is.na(x)])
        x <- as.data.frame(data)
        y <- as.vector(as.data.frame(apply(x, 2, as.numeric)))
        y <- as.vector(unlist(y))
        if (suppressWarnings(length(y[is.na(y)]) - na > 0)) stop("Data contains non-numeric data. Make sure not to use ',' as a decimal operator")
        ## check all data.frames are named
        if ((is.null(names(data)))) {
            pass <- FALSE
        warning("Every data.frame needs to be named with the sample id")
        }
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        if (any(table(col_names) != length(data))) {
            pass <- FALSE
            warning("Every sample needs to have the same number of columns")
        }
        col_names <- names(data[[1]])
        ind_names <- names(data)
        if (any(duplicated(ind_names))) warning("Avoid duplicates in sample names")
        ## Check for presence of any "," instead of a "." as decimal sign
        x <- as.vector(unlist(data))
        na <- length(x[is.na(x)])
        y <- as.vector(as.data.frame(apply(as.data.frame(data), 2, as.numeric)))
        y <- as.vector(unlist(y))
        if (length(y[is.na(y)]) - na > 0) stop("Data contains non-numeric data. Make sure not to use ',' as a decimal operator")
        gc_peak_list <- data
    }

## Do some internal checks if write_output is defined in align_chromatograms
        if (any(names(opt) == "write_output")) {
        if (any(!(opt[["write_output"]] %in% col_names))) {
            pass <- FALSE
            warning("Names in write_output have to be included as a variable in the data!")
        }
            }
    if (any(stringr::str_detect(string = ind_names, pattern = " "))) warning("Avoid whitespaces in Sample Names!")
    if (any(stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Sample Names should only contain Letters, Numbers and '_' ")
    if (any(stringr::str_detect(string = col_names, pattern = " "))) warning("Avoid whitespaces in Variable Names!")
    if (any(stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Variable Names should only contain Letters, Numbers and '_' ")

    if (any(names(opt) == "blank")) {
        if (any(!(opt[["blank"]] %in% ind_names))) {
            pass <- FALSE
        warning("blanks have to refer to samples in the data!")
        }
    }
    if (any(names(opt) == "reference")) {
        if (!is.null(opt[["reference"]]) & any(!(opt[["reference"]] %in% ind_names))) {
    pass <- FALSE
    warning("reference has to be included as a sample in the data!")
}
}
    format_error <- function(x){
        pass <- TRUE
        check_var_count <- function(x) {
            mat <- as.matrix(x)
            L <- length(unique(colSums(!is.na(x))) == 1)
            return(L)
        }
        y <- unlist(lapply(gc_peak_list,check_var_count))
        if (length(which(y != 1)) > 0) {
            out <- names(y[which(y != 1)])
            warning(paste(out,collapse = "; ") ," violate(s) the requirements.",call. = FALSE)
            warning("Every sample needs to have the same number of values for each variable!",call. = FALSE)
            pass <- FALSE
        }
        return(pass)
    }

    ## Checks that every sample has the same number of values per column
    format_pass <- format_error(gc_peak_list)
    if (pass == TRUE & format_pass == TRUE) {
        cat("All checks passed!\nReady for processing with align_chromatograms\n")
    } else {
        cat("Not all checks have been passed. Read warning messages and change data accordingly\n")
    }


    if (plot == TRUE) {
        counter <- function(gc_peak_list){
            number <- lapply(gc_peak_list, function(x){
                ## vectorize the first column
                temp <- x[,1]
                ## Estimate number of peaks
                length(temp[!is.na(temp) & temp > 0])
            } )
            out <- t(as.data.frame((number)))
            out <- reshape2::melt(out)
            out <- out[,c("Var1","value")]
            out <- as.data.frame(out)
            colnames(out) <- c("ID","Peaks")
            return(out)
        }

        out <- counter(gc_peak_list)
        peaks <- as.vector(unlist(out["Peaks"]))
        names(peaks) <- unlist(out["ID"])
        ymax <- max(peaks)

        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = ""))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = ""))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "# Peaks"))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.25))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.25))
        if (!"cex.names" %in% names(mcall)) {
            lab_thresh <- c(20,30,40,50,60,Inf)
            lab_size <- c(1.2,1.1,0.95,0.85,0.75,0.7)
            samples_size <- length(peaks)
            # find the matching size
            temp <- which(lab_thresh > samples_size)
            if (min(temp) == 1) {
                label_size <- lab_size[1]
            } else {
                label_size <- lab_size[min(temp) - 1]
            }
        }
        if (!"cex.names" %in% names(mcall)) arg_list <- append(arg_list,list(cex.names = label_size))

        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "blue"))
        if (!"srt" %in% names(mcall))  arg_list <- append(arg_list,list(srt = 45))
        if (!"las" %in% names(mcall))  arg_list <- append(arg_list,list(las = 2))
        if (!"names.arg" %in% names(mcall)) arg_list <- append(arg_list,list(names.arg = names(peaks)))
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list,list(ylim = c(0,ymax + 5)))

        bars <- do.call(graphics::barplot,args = c(list(height = peaks),arg_list,...))
    }
    if (plot == TRUE) {
        ## Form a data frame with peak numbers
        output <- data.frame(sample = names(peaks), peaks = peaks,row.names = 1:length(peaks))
    } else {
        ## Give a True for fast checks without plotting
        if (pass == TRUE & format_pass == TRUE) {
            output <- "TRUE"
        } else {
            output <- "FALSE"
        }
    }
    ## return
    return.values <- output
}
