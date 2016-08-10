#' Check correct formatting of data input
#'
#'@description
#' Checks conformity between the input format and the requirements of GCalignR. Supported are
#' \code{.txt} files or lists of data frames. See \code{\link{align_chromatograms}} for details.
#'
#'@param data
#'       path to a data file or the name of a list in the Global Environment.
#'
#'@param show_peaks
#'logical, if TRUE the distribution of peak numbers is plotted. Default is FALSE.
#'
#'@param sep
#'The field separator character. Values on each line of the file are separated by this
#'character. The default is tab seperated (sep = '\\t'). See \code{sep} argument in \code{\link[utils]{read.table}} for details.
#'
#'@param ...
#'optional arguments passed to methods, see \code{\link[graphics]{barplot}}. Only used if
#'\code{show_peaks==TRUE}.
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@import magrittr stringr
#'
#' @examples
#' check_input(seal_peaks) ## Checks format
#' check_input(seal_peaks, show_peaks=TRUE) ## Includes a barplot of peaks
#'
#' @export
#'

check_input <- function(data,show_peaks=FALSE, sep = "\t",...) {

    mcall = as.list(match.call())[-1L]
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
        if(any(!(opt[["write_output"]]%in%col_names))) stop("Names in write_output have to be included as a variable in the data!")
    }
    if(any(stringr::str_detect(string = ind_names, pattern = " "))) warning("Avoid whitespaces in Sample Names!")
    if(any(stringr::str_detect(string = ind_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Sample Names should only contain Letters, Numbers and '_' ")
    if(any(stringr::str_detect(string = col_names, pattern = " "))) warning("Avoid whitespaces in Variable Names!")
    if(any(stringr::str_detect(string = col_names, pattern = "[^a-zA-Z\\d\\_]"))) warning("Variable Names should only contain Letters, Numbers and '_' ")

    if(any(names(opt)=="blank")){
        if(any(!(opt[["blank"]]%in%ind_names))) stop("blanks have to refer to samples in the data!")
    }
    if(any(names(opt)=="reference")){
        if(any(!(opt[["reference"]]%in%ind_names))) stop("reference has to be included as a sample in the data!")
    }
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

    if(show_peaks==TRUE){
        counter <- function(gc_peak_list){
            number <- lapply(gc_peak_list, function(x){
                temp <- x[,1] # vectorize the first column
                length(temp[!is.na(temp)]) # number of peaks
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
        if(!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Peak Counts in the raw data\n befire running GCalignR"))
        if(!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = ""))
        if(!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Number of Peaks"))
        if(!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis=1.5))
        if(!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab=1.5))
        if(!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "blue"))
        if(!"srt"%in% names(mcall))  arg_list <- append(arg_list,list(srt = 45))
        if(!"las"%in% names(mcall))  arg_list <- append(arg_list,list(las = 2))
        if(!"names.arg" %in% names(mcall)) arg_list <- append(arg_list,list(names.arg=names(peaks)))
        if(!"ylim" %in% names(mcall)) arg_list <- append(arg_list,list(ylim=c(0,ymax+5)))

        bars <- do.call(graphics::barplot,args=c(list(height=peaks),arg_list,...))
        graphics::text(x=bars,y=peaks+2,labels = as.character(peaks),cex = 0.9)

    }

}


