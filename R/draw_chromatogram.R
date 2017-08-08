#' Visualise peak lists as a chromatogram
#'
#' @description Creates a graphical representation of one or multiple peak lists in the form of simplified chromatograms. Peaks are represented by gaussian distributions centred at the peak retention time.
#'
#' @details Peaks from the peak list are depicted as gaussian distributions. By default, peaks are visualised regardless of their concentration, thereby indicating presence/absence only. If the data is an "GCalign" object that was processed with \code{\link{align_chromatograms}}, chromatograms can be drawn for the dataset prior to alignment (\strong{"pre_aligned"}), after correcting linear drift (\strong{"lin_aligned"}) or after the complete alignment was conducted (\strong{"fully_aligned"}). In the latter case, retention times refer to the mean retention time of a substance scored among samples and do not represent between-sample variation anymore. Depending on the range of retention times and the distance among substances the peak width can be adjusted to enable a better visual separation of peaks by changing the value of parameter \code{width}. Note, aligned peaks (= exactly matching retention time) will overlap completely and only the last sample plotted will be visible. Hence, the number of samples is given on top of each peak but default. The function returns a list containing the ggplot object along with the underlying data frame to allow for maximum control in adapting the plot (see examples section in this document).
#'
#' @param data Currently three formats are supported. Either a \strong{text file} or a \strong{list} containing raw data peak lists or an \strong{GCalign object}. See \code{\link{align_chromatograms}} for details.
#'
#' @inheritParams align_chromatograms
#'
#' @param width numeric value giving the standard deviation of gaussian peaks. Decrease this value to separate overlapping peaks within samples. Default is 0.01.
#'
#' @param step Either \strong{"pre_aligned"}, \strong{"lin_aligned"} or \strong{"fully_aligned"}. Only supported when \code{data} is an object of type "GCalign". By default aligned is used.
#'
#' @param breaks A vector giving the breakpoints between ticks on the x axis.
#'
#' @param rt_limits A vector of length two giving min and max values or retention times to plot.
#'
#' @param samples A character vector of sample names to draw chromatograms of a subset. By default every sample will be drawn.
#'
#' @param show_num Logical indicating whether sample numbers are drawn on top of each peak. By default TRUE.
#'
#' @param show_rt Logical indicating whether peak retention times are drawn on top of each peak. By default FALSE.
#'
#' @param plot Logical indicating if a plot is printed. By default TRUE.
#'
#' @param shape
#' a character determining the shape of peaks. Peaks are approximated as "gaussian" by default. Alternatively, peaks can be visualised as "sticks".
#'
#' @param legend.position
#' See \code{\link[ggplot2]{theme}} for options of legend positions.
#'
#' @return A list containing the data frame created for plotting and the ggplot object. See \code{\link[ggplot2]{ggplot}}.
#'
#' @examples
#' ## load data
#' path <- (system.file("extdata", "simulated_peak_data.txt", package = "GCalignR"))
#' ## run with defaults
#' x <- draw_chromatogram(data = path, rt_col_name = "rt")
#' ## Customise and split samples in panels
#' x <- draw_chromatogram(data = path, rt_col_name = "rt", samples = c("A2","A4"),
#'  plot = FALSE, show_num = FALSE)
#' x[["ggplot"]] + ggplot2::facet_wrap(~ sample, nrow = 2)
#' ## plot without numbers
#' x <- draw_chromatogram(data = path, show_num = FALSE, rt_col_name = "rt")
#'
#' @export
#'
draw_chromatogram <- function(data = NULL, rt_col_name = NULL, width = 0.1, step = NULL, sep = "\t", breaks = NULL, rt_limits = NULL, samples = NULL, show_num = TRUE, show_rt = FALSE, plot = TRUE, shape = c("gaussian","stick"), legend.position = "bottom")  {


# Internal functions
# ##################
rt_min_max <- function(df, rt_col_name) data.frame(min = min(df[[rt_col_name]][df[[rt_col_name]] > 0], na.rm = T), max = max(df[[rt_col_name]], na.rm = T))

conc_max <- function(df, conc_col_name) data.frame(min = min(df[[conc_col_name]][df[[conc_col_name]] > 0]), max = max(df[[conc_col_name]]))

p2c <- function(df = NULL, x = NULL, conc_max = NULL, rt_col_name = NULL, conc_col_name = NULL) {
    temp <- rep(0, length(x))
    if (is.null(conc_max)) for (i in 1:nrow(df)) temp <- temp + dnorm(x, mean = df[[rt_col_name]][i], sd = width)
    if (!is.null(conc_max)) for (i in 1:nrow(df)) temp <- temp + dnorm(x, mean = df[[rt_col_name]][i], 1.1 - df[[conc_col_name]][i]/conc_max)
    return(temp)
}

time_cut <- function(df, rt_col_name, rt_limits) {
    min <- min(rt_limits)
    max <- max(rt_limits)
    df <- subset(df, (df[[rt_col_name]] >= min) & (df[[rt_col_name]] <= max))
return(df)
}

# #################


# some initial checks
# ###################
    if (is.null(data)) stop("Specify 'data'")
    if (is.null(rt_col_name)) stop("Specify 'rt_col_name'")
    if (show_num == TRUE & show_rt == TRUE) {
        warning("Cannot simulataneously annotate peaks with retention time and sample cout. Set show_rt or show_num to FALSE")
        show_rt <- FALSE
    }
shape <- match.arg(shape)


    if (is.character(data)) {
        out <- check_input(data = data, rt_col_name = rt_col_name, sep = sep, plot = F, message = F)
        if (out == FALSE) stop("Data is malformed. See check_input for details")
    } else {
        if (class(data) == "GCalign") {
            if (!(rt_col_name %in% names(data[["aligned"]])))  stop(print(paste(rt_col_name,"is not a valid variable name. Data contains:",paste(names(data[["aligned"]]),collapse = " & "))))
        } else if (class(data) == "list") {
            out <- check_input(data = data, rt_col_name = rt_col_name, sep = sep, plot = F, message = F)
            if (out == FALSE) stop("Data is malformed. See check_input for details")
        }
    }
# ###################

    # data is text file
    if (is.character(data)) {
    # gather data and information
peak_list <- read_peak_list(data, sep, rt_col_name)
    } else if (class(data) == "GCalign") {
        step <- match.arg(step, choices = c("fully_aligned","pre_aligned","lin_aligned"))
        if (step == "pre_aligned") {
            peak_list <- data[["input_list"]]
        } else if (step == "lin_aligned") {
            shift <- data[["Logfile"]][["LinearShift"]]
            peak_list <- data[["input_list"]]
            for (i in 1:length(peak_list)) peak_list[[i]][[rt_col_name]] <- peak_list[[i]][[rt_col_name]] + shift[i,1]
        } else if (step == "fully_aligned") {
            peak_list <- data[["aligned_list"]]
            substance_rts <- data[["aligned"]][[rt_col_name]][["mean_RT"]]
            peak_list <- lapply(X = peak_list, FUN = function(x) {
                x[[rt_col_name]][x[[rt_col_name ]] > 0] <- substance_rts[which(x[[rt_col_name]] > 0)]
                return(x)
        })
        }
    } else if (class(data) == "list") {
        peak_list <- lapply(data, FUN = function(x) {
            if (any(is.na(rowSums(x)))) {
                p <- as.vector(which(is.na(rowSums(x))))
                x <- x[-p,]
            }
            if (any(rowSums(x) == 0)) {
                p <- as.vector(which(rowSums(x) == 0))
                x <- x[-p]
            }
            return(x)
        })
        peak_list <- data
    }

peak_list <- lapply(peak_list, FUN = function(x) {
    if (any(is.na(rowSums(x)))) {
        p <- as.vector(which(is.na(rowSums(x))))
        x <- x[-p,]
    }
    if (any(rowSums(x) == 0)) {
        p <- as.vector(which(rowSums(x) == 0))
        x <- x[-p,]
    }
    return(x)
})
    # remove zeros
peak_list <- lapply(X = peak_list, FUN = function(x) {
    # get empty rows
    r <- which(x[[rt_col_name]] <= 0)
    # remove empty rows when present
    if (length(r) > 0) x[-r,]
    return(x)
} )

    if (!is.null(samples)) peak_list <- peak_list[samples]
samples <- names(peak_list)
if (!is.null(rt_limits)) peak_list <- lapply(peak_list, time_cut, rt_col_name = rt_col_name, rt_limits = rt_limits)
# kick out empty samples
empty <- which(lapply(peak_list, nrow) == 0)
peak_list[empty] <- NULL
if (length(peak_list) == 0) stop("No peaks detected. Are all parameters in the function call plausible?")

temp <- as.vector(unlist(lapply(X = peak_list, FUN = rt_min_max, rt_col_name)))
rt_range <- c(floor(min(temp, na.rm = T)), round(max(temp, na.rm = T)))
# treat corner more nicely
margins <- diff(rt_range)*0.05
#rt_range <- seq(from = rt_range[1] - margins, to = rt_range[2] + margins, length = 10000)
rt_range <- seq(from = 0, to = rt_range[2] + margins, length = 10000)

# calculate chromatograms. For now only looking at retention times
conc_col_name <- NULL
if (!is.null(conc_col_name)) {
    if (!any(colnames(peak_list[[1]]) %in% conc_col_name)) stop(paste0("can not access '",conc_col_name," '. Ensure this is a valid variable name!"))
    conc_max <- max(as.vector(unlist(lapply(X = peak_list, FUN = conc_max, conc_col_name))))
} else {
    conc_max <- NULL
}
x <- rep(rt_range, length(samples))
cat("Computing chromatograms ...\n")

pbapply::pboptions(type = "timer", char = "+", style = 1)
y <- as.vector(unlist(pbapply::pblapply(X = peak_list, FUN = p2c, x = rt_range, rt_col_name, conc_col_name = conc_col_name, conc_max = conc_max)))
# y <- as.vector(unlist(lapply(X = peak_list, FUN = p2c, x = rt_range, rt_col_name, conc_col_name = conc_col_name, conc_max = conc_max)))
sample <- rep(samples, each = length(rt_range))
if (isTRUE(show_num)) {
    y2 = rep(0, length(y))
    n <- y2
    df <- data.frame(x,y,sample, y2 = y2, n = n)
} else {
    df <- data.frame(x,y,sample)
}


# plot
if (is.null(breaks)) breaks <- seq(from = min(rt_range), to = max(rt_range), by = 5)

## detect peaks for all chromatograms
peaks <- find_peaks(df)
min_x <- ifelse(min(peaks[["x"]]) - margins > 0,min(peaks[["x"]]) - margins,0)
df <- subset(df, x > min_x)

if (shape == "gaussian") {
chroma <- ggplot(data = df, aes(x,y, col = sample)) + geom_line(size = 1)
} else if (shape == "stick") {
    chroma <- ggplot(data = peaks, aes(x, y, col = sample)) + geom_point() + geom_linerange(data = peaks, aes(x = x, ymin = 0, ymax = y), linetype = "dashed")
}

chroma <- chroma + theme_classic() + xlab("Retention time") + ylab("") + scale_x_continuous(breaks = breaks, expand = c(0,0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = legend.position) + guides(col = guide_legend(ncol = 10, title = NULL))

if (isTRUE(show_num)) {
## count number of substances per peak
pn <- data.frame(x = as.numeric(names(summary(as.factor(peaks[["x"]]), maxsum = length(peaks[["x"]])))), n = as.vector(summary(as.factor(peaks[["x"]]), maxsum = length(peaks[["x"]]))))
pn[["y"]] <- unlist(lapply(pn[["x"]], FUN = function(x) {
    max(peaks[["y"]][round(peaks[["x"]],4) == round(x,4)])
    }))

for (i in 1:nrow(pn)) {
    tx <- which(round(df[["x"]],4) == round(pn[["x"]][i],4))[1]
    df[["y2"]][tx] <- pn[["y"]][i]
    df[["n"]][tx] <- pn[["n"]][i]
}

chroma <- chroma +
     geom_segment(data = subset(df, n > 0), aes(x = x, xend = x, y = y2, yend = 0), colour = "black", size = 0.5, linetype = "dotted") +
    annotate("text", x = pn[["x"]] , y = pn[["y"]] + 0.5, label = as.character(pn[["n"]]))
}

if (isTRUE(show_rt)) {
    peaks <- find_peaks(df)
    chroma <- chroma + geom_linerange(data = peaks, aes(x = x, ymin = y, ymax = y + 0.1), linetype = "solid", col = "black") + annotate("text", x = peaks[["x"]], y = peaks[["y"]] + 0.25, label = as.character(round(peaks[["x"]],2)), angle = 90)
}

cat("\nDrawing chromatograms ...")
if (plot == TRUE) print(chroma)

cat("\nDone\n")

## output
###########
output <- list(ggplot = chroma, df = df)
###########

}#end draw_chroma
