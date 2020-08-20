#' Estimate the observed space between peaks within chromatograms
#'
#'@description
#' The parameter \code{min_diff_peak2peak} is a major determinant in the alignment of a dataset with \code{\link{align_chromatograms}}. This function helps to infer a suitable value based on the input data. The underlying assumption here is that distinct peaks within a separated by a larger gap than homologous peaks across samples. Tightly spaced peaks within a sample will appear on the left side of the plotted distribution and can indicate the presence of split peaks in the data.
#'
#' @inheritParams check_input
#'
#' @inheritParams align_chromatograms
#'
#' @param quantile_range
#' A numeric vector of length two that allows to subset an arbitrary interquartile range.
#'
#' @param quantiles
#' A numeric vector. Specified quantiles are calculated from the distribution.
#'
#' @param by_sample
#' A logical that allows to calculate peak interspaces individually for each sample. By default all samples are combined to give the global distribution of next-peak differences in retention times. When \code{by_sample = TRUE}, a series of plots (one for each sample) is created and a keystroke is required to proceed.
#'
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
peak_interspace <- function(data,rt_col_name = NULL, sep = "\t", quantiles = NULL, quantile_range = c(0,1), by_sample = FALSE) {
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
        gc_peak_list <- data

        ## remove rows only containing NA.
        gc_peak_list <- lapply(gc_peak_list, function(gc_data) {
            gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
            ## Remove empty rows
            gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
            as.data.frame(gc_data)
            })
        ## coerce to numeric
        gc_peak_list <- lapply(gc_peak_list, function(x) {
            x[[rt_col_name]] <- as.numeric(x[[rt_col_name]])
            return(x)
        })
    }



    # Extract retention times for all samples and list them
    fx <- function(df,rt_col_name) out <- df[[rt_col_name]]

    rt_list <- lapply(gc_peak_list,FUN = fx,rt_col_name = rt_col_name)
    # Calculate peak interspaces for each sample
    spaces <- round(unlist(lapply(rt_list,FUN = function(x) diff(x[x > 0 & !is.na(x)]))),2)
    spaces_table <- table(as.factor(spaces))
    breaks <- as.numeric(as.character(names(spaces_table)))
# Prepare barplot/histogram
    bar_data <- spaces_table[min(which(breaks >= stats::quantile(spaces,quantile_range[1]))):min(which(breaks >= stats::quantile(spaces,quantile_range[2])))]

    if (by_sample == FALSE) {
        graphics::plot(bar_data/sum(bar_data),xpd = T,col = "darkblue",xlab = "Distance between neighbouring peaks\nwithin samples",ylab = "Proportion")
    } else if (by_sample == TRUE) {
        dat <- lapply(gc_peak_list, FUN = fx, rt_col_name = rt_col_name)
        names <- names(dat)
        dat_spaces <- lapply(dat, function(x) round(diff(x[x > 0 & !is.na(x)]), 2))
        dat_table <- lapply(dat_spaces, function(x) table(as.factor(x)))
        escape <- "start"
        while (!(escape %in% c("stop","Stop","escape","ESC","Esc","Escpape"))) {
        for (i in 1:length(dat_table)) {
            graphics::plot(dat_table[[i]]/sum(dat_table[[i]]),xpd = T,col = "darkblue",xlab = "Distance between neighbouring peaks\nwithin samples",ylab = "Proportion", main = paste0(names[i], " (", i, "/", length(dat_table),")"))
          escape <-  readline("Press enter to proceed or Esc to interrupt: ")
          cat(escape)
        }
        }
    }

    spaces_subset <- spaces[spaces >= stats::quantile(spaces,quantile_range[1]) & spaces <= stats::quantile(spaces,quantile_range[2])]


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
