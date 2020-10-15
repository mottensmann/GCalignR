#' Check input prior to processing in GCalignR
#'
#'@description
#' Checks input files for common formatting problems.
#'
#' @details
#' Sample names should contain just letters, numbers and underscores and no whitespaces.
#' Each sample has to contain the same number of columns, one of which is the retention
#' time and the others are arbitrary variables in consistent order across samples. Retention times are expected to be numeric, i.e. they are only allowed to contain numbers from 0-9 and "." as the only decimal character. Have a look at the vignettes for examples.
#'
#'@param plot
#' Boolean specifying if the distribution of peak numbers is plotted.
#'
#'@inheritParams align_chromatograms
#'
#'@param ...
#'optional arguments passed to methods, see \code{\link[graphics]{barplot}}.
#'
#'@param message
#'Boolean determining if passing all checks is indicated by a message.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @examples
#' ## gc-data
#' data("peak_data")
#' ## Checks format
#' check_input(peak_data)
#' ## Includes a barplot of peak numbers in the raw data
#' check_input(peak_data, plot = TRUE)
#'
#' @export
#'
check_input <- function(data,plot = FALSE, sep = "\t", message = TRUE, ...) {

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
        ind_names <- readr::read_lines(data, n_max = 1)
        ind_names <- unlist(stringr::str_split(string = ind_names,pattern = sep))
        ind_names <- ind_names[ind_names != ""]
        ## Variable Names
        col_names <- readr::read_lines(data, n_max = 1, skip = 1)
        col_names <- unlist(stringr::str_split(string = col_names,pattern = sep))
        col_names <- col_names[col_names != ""]
        col_names <- stringr::str_trim(col_names)
        # validate retention time name
        if ("rt_col_name" %in% names(opt)) {
            rt_col_name <- opt[["rt_col_name"]]
            if (!(rt_col_name %in% col_names)) {
                stop(print(paste(rt_col_name,"is not a valid variable name. Data contains:",paste(col_names,collapse = " & "))))
                pass <- FALSE
            }
            # check conc_col_name
            if ("conc_col_name" %in% names(opt)) {
                conc_col_name <- opt[["conc_col_name"]]
                if (!(conc_col_name %in% col_names)) {
                    stop(print(paste(conc_col_name,"is not a valid variable name. Data contains:",paste(col_names, collapse = " & "))))
                    pass <- FALSE
                }
            }
        }
        ind_names <- stringr::str_trim(ind_names)
        ## Get Peak Data
        gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F, fill = T)


        if (!"rt_col_name" %in% names(opt)) {
            col_class <- as.factor(unlist(lapply(lapply(X = 1:ncol(gc_data), function(x) as.vector(gc_data[,x])), class)))
            if (any(!(col_class %in% c("numeric", "integer")))) {
                message(
                    paste0("Retention times need to be numeric. Column(s)",
                           paste(as.character(which(!col_class %in% c("numeric", "integer"))), collapse = "; "),
                           " violate the requirements." ))
            }
        }


        ## Remove just NA-rows
        gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
        ## Remove empty rows
        gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]

        gc_data <-  as.data.frame(gc_data)

        ### Check input for completeness ###

        if (!((ncol(gc_data) / length(col_names)) %% 1) == 0) {
            pass <- FALSE
            stop("Number of data columns is not a multiple of the column names provided")
        }
        if (!((ncol(gc_data) / length(col_names))  == length(ind_names))) {
            pass <- FALSE
            stop(paste0("Number of sample names provided does not fit to the number of columns in the data."),
                 "\n\tExpected # = ", ncol(gc_data) / length(col_names), " ;Provided # = ", length(ind_names))
        }
        if (any(duplicated(ind_names))) {
            warning(paste0("Duplicated sample names are not allowed.\nChange sample(s):\n",as.character(ind_names[duplicated(ind_names)])))
            pass <- FALSE
        }
        ## convert to list
        gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)


        ## If data is a list of data.frames
    } else if (is.list(data)) {
        ## check that every element in the list is a data.frame
        if (!(all(unlist(lapply(data, is.data.frame))))) {
            pass <- FALSE
            warning("Every Sample has to be a data.frame")
        } else if (all(unlist(lapply(data, tibble::is_tibble)))) {
            message("Detected tibbles, coerce to data frame")
            data <- lapply(data, as.data.frame)
        }

        # check class of columns
        min_n <- min(as.vector(unlist(lapply(data, nrow))))
        data2 <- lapply(data, function(x) x[1:min_n,])
        temp <- do.call("cbind", data2)
        if (!(any(apply(temp, 2, class) %in% c("numeric","integer")))) {
            # na_1 <- length(which(is.na(temp)))
            # data <- lapply(data, function(x) as.data.frame(apply(x, 2, as.numeric)))
            # data2 <- lapply(data, function(x) x[1:min_n,])
            # na_2 <- length(which(is.na(do.call("cbind", data2))))
            message("Non-numeric variables. Make sure all retention time values are numeric\n")
            # if (na_2 > na_1) {
            #     pass <- FALSE
            #     warning("NAs introduced by coercion\n")
            # }
        }

        ## Adjust sizes of data frames, each has the same number of rows
        data <- lapply(data,matrix_append,gc_peak_list = data, val = "NA")
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
        # Validate retention time variable
        if ("rt_col_name" %in% names(opt)) {
            rt_col_name <- opt[["rt_col_name"]]
            if (!(rt_col_name %in% col_names)) stop(print(paste(rt_col_name,"is not a valid variable name. Data contains:",paste(col_names,collapse = " & "))))
        }
        ind_names <- names(data)
        if (any(duplicated(ind_names))) warning("Avoid duplicates in sample names\n")
        gc_peak_list <- data
    }# end of checking on list

    ## check that each peak consistently is ordered by increasing
    ## peak retention times
    if ("rt_col_name" %in% names(opt)) {
        ordered.input <- sapply(gc_peak_list, function(x) {
            any(diff(order(x[[rt_col_name]])) != 1)
        })
        if (any(ordered.input == TRUE)) {
            warning("At least one peak list contains unordered retention times.\n Retention times will be reordered by increasing retention times, (i.e. from low to high values)\n")
            gc_peak_list <- lapply(gc_peak_list, function(x) {
                ## sort by increasing retention times
                x[order(x[[rt_col_name]], decreasing = F),]
            })
        }
    }


    # Validate retention times
    if ("rt_col_name" %in% names(opt)) {
        df <- unlist(lapply(gc_peak_list,FUN = function(x,rt_col_name) x[[rt_col_name]],rt_col_name))
        if (!is.numeric(df)) stop("Not all retention times are numeric! Make sure to use the correct decimal operator '.' instead of ','")
    }

    if ("conc_col_name" %in% names(opt)) {
        df <- unlist(lapply(gc_peak_list,FUN = function(x,rt_col_name) x[[conc_col_name]],conc_col_name))
        if (!is.numeric(df)) stop("Not all concentration values are numeric! Make sure to use the correct decimal operator '.' instead of ','")
    }


    ## Do some internal checks if write_output is defined in align_chromatograms
    if (any(names(opt) == "write_output")) {
        if (any(!(opt[["write_output"]] %in% col_names))) {
            pass <- FALSE
            warning("Names in write_output have to be included as a variable in the data!")
        }
    }
    ## Lines checking for whitespaces are redundant!

    # if (any(stringr::str_detect(string = ind_names, pattern = " "))) warning("Avoid whitespaces in Sample Names!")

    # Check for proper individual names
    if (any(stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]"))) {
        warning("Avoid whitespaces in sample names! Additionally they should only contain Letters, Numbers and '_' ","\n",paste(ind_names[stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]")],collapse = "; ")," violate(s) these requirements.")
    }
    #if (any(stringr::str_detect(string = col_names, pattern = " "))) warning("Avoid whitespaces in Variable Names!")

    # check for proper variable naming
    if (any(stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]"))) {
        warning("Avoid whitespaces in variable names! Additionally they should only contain Letters, Numbers and '_' ","\n",paste(col_names[stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]")],collapse = "; ")," violate(s) these requirements.")
    }

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
            warning("Every sample needs to have the same number of values for each variable!","\n", paste(out,collapse = "; ")," violate(s) the requirements.",call. = FALSE)
            #warning(paste(out,collapse = "; ") ," violate(s) the requirements.",call. = FALSE)
            pass <- FALSE
        }
        return(pass)
    }

    ## Checks that every sample has the same number of values per column
    format_pass <- format_error(gc_peak_list)
    if (pass == TRUE & format_pass == TRUE) {
        if (message == TRUE) cat("All checks passed!\n\n")
    } else {
        cat("\nNot all checks have been passed. Read warning messages below and change accordingly to proceed\n\n")
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
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Number of Peaks"))
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
