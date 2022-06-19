#' Remove singletons
#'
#' @description
#' Identifies and removes singletons (i.e. peaks that are unique for one sample) from the aligned dataset.
#'
#' @inheritParams remove_blanks
#'
#' @return a list of data frames for each individual.
#'
#' @examples
#' data("peak_data")
#' ## subset for faster processing
#' data <- lapply(peak_data[1:5], function(x) x[20:35,])
#' x <- align_chromatograms(data, rt_col_name = "time")
#' out <- remove_singletons(data = x)
#'
#' @author Meinolf Ottensmann (meinolf.ottensmann@@web.de) & Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
remove_singletons <- function(data) {
    if (inherits(data, "GCalign")) {
        rt_col_name <-  data[["Logfile"]][["Call"]][["rt_col_name"]]
        data <- data[["aligned_list"]]
    } else if (is.list(data)) {
        # does not matter which column to take
        rt_col_name <- names(data[[1]])[1]
    } else {
        stop("data is not of class GCalign or a list of data frames")
    }

    rt_mat <- do.call(cbind, lapply(data, function(x) x[[rt_col_name]]))
    indices <- which(rowSums(rt_mat > 0) == 1)
    if (length(indices) > 0) {
        data <- lapply(data, function(x) x[-indices, ])
    }
    cat(paste(as.character(length(indices)),'singletons were removed\n'))
    return(data)
}
