#' Subtraction of blank readings from sample readings
#'
#' @description For each substance that is present in blanks, samples are corrected by subtraction of the respective quantity. If more than one sample is submitted, abundances are averaged. This procedure is sensitive to differences in the total concentration of samples and should be applied to samples where the preparation yields comparable concentrations for each sample.
#'
#' @details Substances that are present in one or more blanks are identified in the aligned dataset, then the mean abundance is calculated for the blanks and the corresponding value is subtracted from each sample. If the control contains higher concentration (i.e. blank substation creates negative abundances) warnings will be shown and the respective value will be set to zero
#'
#' @param blanks
#' Character vector of names of negative controls.
#'
#' @param input
#' A GCalign Object
#'
#' @param conc_col_name
#' If the input is a GCalign object the variable containing the abundance values needs to be specified.
#'
#' @keywords beta
#'
#' @examples
#' ## Not run
#' #out <- blank_substraction(aligned_peak_data, blanks = "M2", conc_col_name = "area")
#'
#'
blank_substraction <- function(input = NULL, blanks = NULL, conc_col_name = NULL) {

    # Internal functions
    # -------------------------------------------
    get_contaminants <- function(blanks, data) {
        # find peak indices
        fx <- function(x, y) which(y[[x]] > 0)
        out <- lapply(X = as.list(blanks),FUN = fx, y = data)
        return(sort(unique(unlist(out))))
    }

    name_in_data <- function(name, input) {
        if (!(name %in% names(input))) (paste(as.character(name), "is not available as a sample name"))
    }
    calc_conc <- function(row, blanks, input) {
        df <- as.data.frame(input)
        return(mean_no_zero(df[row,unlist(blanks)]))
    }

    mean_no_zero <- function(x) mean(x[x != 0])

    blank_substract <- function(input = input, substances = unlist(rows), conc = conc) {

        fx <- function(x) {
            temp <- x[substances] - unlist(conc)
            if (any(temp < 0)) warning("Negative values introduced by blank substraction")
            temp[temp < 0] <- 0
            x[substances] <- temp
            return(x)
        }
        return(lapply(input,fx))
    }
    # -------------------------------------------

# Prepare input data and do some checks
if (is.null(input)) stop("No input was defined")
if (is.null(blanks)) stop("Define name(s) of blanks")

    # read data and prepare a list
if (inherits(input, "GCalign")) {
    if (is.null(conc_col_name)) stop("Define the name of a data frame")
    if (conc_col_name %in%  names(input[["aligned"]])) {
    input2 <- input[["aligned"]][[conc_col_name]]
    } else {
        stop("conc_col_name is not a valid name of a data frame")
        }
    input2 <- as.list(input2)[-1]
}
    # ensure blanks are samples of data
out <- unlist(lapply(as.list(blanks), name_in_data, input2))
if (length(out) > 1) stop(cat(paste(out,collapse = "\n")),"Issue with blank names")
blank_names <- blanks
blanks <- as.list(blanks)

    # determine contaminations
rows <- as.list(get_contaminants(blanks, data = input2))

    # calculate mean concentration of controls
conc <- lapply(rows, calc_conc, blanks = blanks, input = input2)

    # subtract background concentration from samples
out <- blank_substract(input = input2, substances = unlist(rows),conc = conc)

## put back into GCalign object
for (i in 1:length(out)) {
    input[["aligned"]][[conc_col_name]][,i + 1] <- out[[i]]
    input[["aligned_list"]][[i]][[conc_col_name]] <- out[[i]]
}

## Remove blanks
index <- which(names(input[["aligned"]][[conc_col_name]]) %in% blank_names)
for (i in 1:length(input[["aligned"]])) {
    input[["aligned"]][[i]] <- input[["aligned"]][[i]][,-index]

}
input[["aligned_list"]][blank_names] <- NULL


return(input)
}
