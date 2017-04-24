#' Estimate the observed space between peaks within chromatograms
#'
#'@description
#' The parameter \code{min_diff_peak2peak} is a major determinant in the alignment of a dataset with \code{\link{align_chromatograms}}.
#' This function allows to infer a suitable value based on the input data.
#'
#' @inheritParams check_input
#' @inheritParams align_chromatograms
#' @param quantile_range
#' A numeric vector of length two specifying arbitrary interquartile ranges visualised in a barplot. By default the whole range is plotted \code{quantile_range = c(0,1)}.
#' @param quantiles
#' A numeric vector. Specified quantiles are calculated from the distribution.
#' @return List containing summary statistics of the peak interspace distribution
#' @import stringr
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @examples
#' ## plotting with defaults
#' peak_interspace(data = peak_data, rt_col_name = "time")
#' ## plotting up to the 0.95 quantile
#' peak_interspace(data = peak_data,rt_col_name = "time",quantile_range = c(0,0.95))
#' ## return the 0.1 quantile
#' peak_interspace(data = peak_data,rt_col_name = "time", quantiles = 0.1)
#'
#' @export
#'
peak_interspace <- function(data,rt_col_name = NULL, sep = "\t", quantiles = NULL, quantile_range = c(0,1)) {
# Checks
# ######
    if (is.null(rt_col_name)) stop("Specify rt_col_name")
    if (is.character(data)) {
    a <- utils::capture.output(check_input(data = data, sep = sep,rt_col_name = rt_col_name,plot = F))
    if (missing(a)) check_input(data = data, sep = sep,rt_col_name = rt_col_name,plot = F)
    } else {
    a <- utils::capture.output(check_input(data = data,rt_col_name = rt_col_name,plot = F))
    if (missing(a)) check_input(data = data,rt_col_name = rt_col_name,plot = F)
    }
# 2. Load Data
# ############
    if (is.character(data)) {
        ind_names <- readr::read_lines(data, n_max = 1)
        ind_names <- unlist(stringr::str_split(string = ind_names,pattern = sep))
        ind_names <- ind_names[ind_names != ""]
        ind_names <- stringr::str_trim(ind_names)
        col_names <- readr::read_lines(data, n_max = 1, skip = 1)
        col_names <- unlist(stringr::str_split(string = col_names,pattern = sep))
        col_names <- col_names[col_names != ""]
        col_names <- stringr::str_trim(col_names)

        gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F)
        gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
        gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
        gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))
        gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)

    } else if (is.list(data)) {
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        col_names <- names(data[[1]])
        ind_names <- names(data)
        gc_peak_list <- data }
    # Extract retention times for all samples and list them
    fx <- function(df,rt_col_name) out <- df[[rt_col_name]]
    rt_list <- lapply(gc_peak_list,FUN = fx,rt_col_name = rt_col_name)
    # Calculate peak interspaces for each sample
    spaces <- round(unlist(lapply(rt_list,FUN = function(x) diff(x[x > 0 & !is.na(x)]))),2)
# Put interspaces in bins
    spaces_table <- table(as.factor(spaces))
    breaks <- as.numeric(as.character(names(spaces_table)))
# Prepare barplot/histogram
    bar_data <- spaces_table[min(which(breaks >= stats::quantile(spaces,quantile_range[1]))):min(which(breaks >= stats::quantile(spaces,quantile_range[2])))]
    mids <- graphics::barplot(bar_data/sum(bar_data),las = 2,xpd = T,col = "darkblue",xlab = "Peak interspace",ylab = "Proportion of peaks")
    spaces_subset <- spaces[spaces >= stats::quantile(spaces,quantile_range[1]) & spaces <= stats::quantile(spaces,quantile_range[2]) ]
    # Generate a list to capture output
    output <- list(Summary = summary(spaces))
    if (!is.null(quantiles)) {
        x <- as.list(quantiles)
        fq <- function(quantile,sample) {
            stats::quantile(sample,quantile)
            }
        y <- lapply(x,fq,sample = spaces)
        fn <- function(x) names(x)
        names <- unlist(lapply(y,FUN = fn))
        y <- unlist(y)
        names(y) <- names
        output <- append(output,list(Quantiles = y))
    }
    return(output)
} # end function
