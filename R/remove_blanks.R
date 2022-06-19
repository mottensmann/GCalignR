#' Remove peaks present in negative control samples
#'
#' @description
#' Removes peaks that are present in blanks (i.e. negative control samples) to eliminate contaminations in the aligned data. Afterwards, blanks are deleted itself. This function is only applicable when blanks were not discarded during a previous alignment using \code{\link{align_chromatograms}}.
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
#'@author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
remove_blanks <- function(data, blanks) {
    if (inherits(data, "GCalign")) {
        rt_col_name <-  data[["Logfile"]][["Call"]][["rt_col_name"]]
        data <- data[["aligned_list"]]
    } else if (is.list(data)) {
        # does not matter which column to take
        rt_col_name <- names(data[[1]])[1]
    } else {
            stop("data is not of class GCalign or a list of data frames")
    }

    del_peaks <- function(blanks, data) {
        # get indices

        delete <- sort(unique(unlist(lapply(blanks, function(fx) {
            which(data[[fx]][[rt_col_name]] > 0)
        }))))

        # remove blanks
        data[blanks] <- NULL

        # remove peaks from samples
        if (length(delete) > 0) data <- lapply(data, function(x) x[-delete,])
        return(data)
    }
    data <- del_peaks(blanks = blanks, data = data)
return(data)
    }

