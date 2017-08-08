#' Remove contaminations
#'
#' @description Removes peaks that are present in blanks (i.e. negative control samples) to elliminate contaminations in the aligned data. Afterwards, blanks are deleted itself. This function is only applicable when blanks were not discarded during a previous alignment using \code{\link{align_chromatograms}}.
#'
#' @inheritParams align_chromatograms
#'
#' @param data
#' An object of class "GCalign". See \code{\link{align_chromatograms}} for details. Alternatively, a list of data frames. Whereby each data frame contains the peak list for an individual sample.
#'
#' @return a list of data frames for each individual.
#'
#' @examples
#' data("peak_data")
#' ## subset for faster processing
#' data <- lapply(peak_data[1:5], function(x) x[20:35,])
#' x <- align_chromatograms(data, rt_col_name = "time")
#' out <- remove_blanks(data = x, blanks = c("C2","C3"))
#' ## number of deleted peaks
#' nrow(x[["aligned_list"]][["M2"]]) - nrow(out[["M2"]])
#'
#' @export
#'
remove_blanks <- function(data, blanks) {
    if (class(data) == "GCalign") {
        rt_col_name <-  data[["Logfile"]][["Call"]][["rt_col_name"]]
        data <- data[["aligned_list"]]
    } else if (is.list(data)) {
        # does not matter which column to take
        rt_col_name <- names(data[[1]])[1]
    } else {
            stop("data is not of class GCalign or a list of data frames")
    }
    del_peaks <- function(x, data) {
        # get indices
        del_substances <- which(data[[x]][[rt_col_name]] > 0)
        # remove respective rows
        chroma_out <- lapply(data, function(x) x[-del_substances,])
        # remove blanks
        chroma_out[x] <- NULL
        return(chroma_out)
    }
    for (i in blanks) data <- del_peaks(x = i, data = data)
return(data)
    }

