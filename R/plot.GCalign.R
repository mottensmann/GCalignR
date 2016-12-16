#' Plot Diagonstics for an Gcalign Object
#'
#' @description
#' Four plots are currently available: One plot visualises the distribution of linear shifts
#' that were applied in order to align chromatograms to a reference before aligning individual peaks.
#' A second plot illustrates the remaining variation of retention times on the level of individual
#' peaks by plotting the distribution of retention time ranges. The third plots shows a distribution
#' of peak numbers after aligning the chromatograms. A fourth plot illustrates the amount of peak sharing among chromatograms in a histogram.
#'
#' @examples
#' ## All plots are shown by default
#' plot(aligned_peak_data)
#'
#' ## Distribution of Peaks
#' plot(aligned_peak_data,which_plot="peak_numbers")
#'
#' @param x \code{GCalign} object, created by \code{\link{align_chromatograms}}
#'
#' @param which_plot
#' character string indicating which plot is returned. Available are
#' \strong{"shifts"} a histogram of linear adjustments undertaken in aligning chromatograms,
#' \strong{"variation"} a histogram summarising the range of retention times for every peak defined
#' by the difference between minimum and maximum retention times respectively. The third option
#' is \strong{"peak_numbers"} plotting a barchart of the number of sharings per sample. Additionally \strong{"peaks_shared"} produces a histogram of the proportion with which sharings are shared among samples. This means for every sharing the proportion of samples containing the respective peak is estimated. By default all plots are returned as subplots of one figure.
#'
#' @param ...
#' optional arguments passed on to methods. See
#' \code{\link[graphics]{plot}}, \code{\link[graphics]{hist}} and \code{\link[graphics]{barplot}}.
#' Please Note that optional arguments are currently not passed on when plotting all figures.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @return
#' Depending on the value of \code{which_plot} a data frame containing the data source of the respective plot is returned. If \code{which_plot = "all"} no output is returned.
#'
#' @export
#'
plot.GCalign <- function(x,which_plot = c("all","shifts","variation","peak_numbers","peaks_shared"), ...) {

    ## initialising by picking arguments from call to plot.GCalign
    mcall = as.list(match.call())[-1L]
    ## allow partial matching
    which_plot <- match.arg(which_plot)

### Define internal functions
# -------------------------------------------------------------------
    hist_linshift <- function(object,mcall,...){
        ## Get the applied shifts
         xl <- object[["Logfile"]][["LinearShift"]]["shift"]
         ## steps of linear shifts and their frequency
        df <- as.vector(unlist(as.vector(xl)))
        xmax <- object[["Logfile"]][["Call"]][["max_linear_shift"]]
        xmin <- -xmax

        ## check for optional arguments in the function call, take defaults, if missing
        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Linear Shifts"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Shift size"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Frequency [%]"))
        if (!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(from = xmin,to = xmax + 0.01, by = 0.01)))
        if (!("freq") %in% names(mcall) || !("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = FALSE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.5))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.5))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#1b9e77"))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"xaxt" %in% names(mcall)) arg_list <- append(arg_list, list(xaxt = "n"))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))
        # helper to find a good ylim
        p <- as.vector(as.numeric(summary(as.factor(df))))
        p <- max(p)/sum(p)*100
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))

        x <- do.call(graphics::hist,args = c(list(x = df),arg_list,list(...)))
        graphics::axis(side = 1, at = x[["mids"]], labels = seq(xmin, xmax, 0.01))
        return(df)
    }#end hist_linshift

    hist_peakvar <- function(x, mcall, ...) {
        ## Estimate the range of retention times per sharing, they should be no overlapp
    MinMax <- function(x) {
        temp <- matrix(NA,1,2)
        colnames(temp) <- c("range","row")
        data <- temp[0,]
        for (i in 1:ncol(x)) {
            data <- rbind(data,cbind(abs(diff(range(x[,i][x[,i] > 0],na.rm = TRUE))),i))
        }
        return(as.data.frame(data))
    }
    ## Range of RTs aligned
        df <- MinMax(x[["heatmap_input"]][["aligned_rts"]][,-1])
        ## Formatting
        df <- round(unlist(df["range"]),digits = 2)
        xmax <- max(df)
        xmin <- min(df)


        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Variation in retention times\n(Peaks)"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Range [max - min]"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Frequency [%]"))
        if (!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(xmin,xmax + 0.01,by = 0.01)))
        if (!("freq") %in% names(mcall) || !("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = FALSE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.5))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.5))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#d95f02"))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"xaxt" %in% names(mcall)) arg_list <- append(arg_list, list(xaxt = "n"))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))
        # helper to find a good ylim
        p <- as.vector(as.numeric(summary(as.factor(df))))
        p <- max(p)/sum(p)*100
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))
        x <- do.call(graphics::hist,args = c(list(x = df),arg_list,list(...)))
        graphics::axis(side = 1, at = x[["mids"]], labels = seq(xmin, xmax, 0.01))
        return(df)
    }



    bar_peakdistr <- function(x,mcall,...) {
        rt_var_name <- x[["Logfile"]][["Input"]][["Retention_Time"]]
        conc_var_name <- x[["Logfile"]][["Input"]][["Concentration"]]
        ## Peaks of All Samples
        data <- (x[["aligned"]][[rt_var_name]])
        ## get rid of mean retention time column
        data <- data[,2:ncol(data)]

        peak_df <- matrix(NA,ncol = 2,nrow = length(data))
        peak_df[,1] <- names(data)
        peak_df[,2] <- unlist(lapply(1:ncol(data), function(y) temp <- length(data[,y][data[,y] > 0])))
        peak_df <- data.frame(peak_df)
        names(peak_df) <- c("ID","Peaks")
        peak_df[["Peaks"]] <- as.numeric(as.character(peak_df[["Peaks"]]))

        peaks <- peak_df[["Peaks"]]
        names(peaks) <- peak_df[["ID"]]

        ymax <- max(peaks)
        ymin <- min(peaks)

        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Number of peaks after alignment"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = ""))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Peaks"))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.5))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.5))
        if (!"cex.names" %in% names(mcall)) {
             lab_thresh <- c(20,30,40,50,60,Inf)
             lab_size <- c(1.2,1.1,0.95,0.85,0.75,0.7)
             samples_size <- ncol(data)
             # find the matching size
             temp <- which(lab_thresh > samples_size)
             label_size <- lab_size[min(temp) - 1]
        }
        if (!"cex.names" %in% names(mcall)) arg_list <- append(arg_list,list(cex.names = label_size))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#7570b3"))
        if (!"srt" %in% names(mcall))  arg_list <- append(arg_list,list(srt = 45))
        if (!"las" %in% names(mcall))  arg_list <- append(arg_list,list(las = 2))
        if (!"names.arg" %in% names(mcall)) arg_list <- append(arg_list,list(names.arg = ""))
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list,list(ylim = c(0,ymax + 5)))

        bars <- do.call(graphics::barplot,args = c(list(height = peaks),arg_list,list(...)))
        # graphics::text(x = bars,y = peaks + 2,labels = as.character(peaks),cex = 0.9)
        lablist <- as.vector(names(peaks))
        if (!"names.arg" %in% names(mcall)) {
        graphics::text(bars, graphics::par("usr")[1], labels = lablist, srt = 90, pos = 1, xpd = TRUE, cex = label_size)
        }
        return(peaks)
    }
    hist_shared_peaks <- function(x,mcall,...){
        ### check that this is up to date!
        df <- x[["heatmap_input"]][["aligned_rts"]][-1]
        peaks_shared <- data.frame(time = round(as.numeric(names(df)),3),prop = unlist(lapply(1:ncol(df), function(col) {
            ## all entries of a column
            sub <- df[,col]
            # number of samples having a peak
            N <- length(sub[sub > 0])
            # proportion
            N/length(sub)*100
        })),n = unlist(lapply(1:ncol(df), function(col) {
            ## all entries of a column
            sub <- df[,col]
            # number of samples having a peak
            N <- length(sub[sub > 0])
           })))

        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Sharing of substances"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Frequency of samples [%]"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Number of Substances"))
        if (!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(0,101,by = 1)))
        if (!("freq") %in% names(mcall) || !("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = TRUE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.5))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.5))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#e7298a"))
        # helper to find a good ylim
        p <- max(as.vector(summary(as.factor(round(peaks_shared$prop)))))
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"xaxt" %in% names(mcall)) arg_list <- append(arg_list, list(xaxt = "n"))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))

        x <- do.call(graphics::hist,args = c(list(x = peaks_shared[["prop"]]),arg_list,list(...)))
        graphics::axis(side = 1, at = x[["mids"]][seq(1,101,2)], labels = seq(0, 100, 2))
        return(peaks_shared)
    }
# --------------------------------------------------------------------

if (which_plot == "shifts") {
out <- hist_linshift(object = x,mcall = mcall,...)
out <- data.frame(shifts = out)
out
} else if (which_plot == "variation") {
out <- hist_peakvar(x = x,mcall = mcall,...)
out <- data.frame(range = as.vector(out), index_sharing = 1:length(as.vector(out)))
out
} else if (which_plot == "peak_numbers") {
out <- bar_peakdistr(x = x,mcall = mcall,...)
out <- data.frame(sample = names(out), sharings = out,row.names = 1:length(out))
out
} else if (which_plot == "all") {
    graphics::layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
    bar_peakdistr(x = x,mcall = mcall)
    hist_linshift(object = x,mcall = mcall)
    hist_peakvar(x = x,mcall = mcall)
    hist_shared_peaks(x = x, mcall = mcall,...)
    ## back to normal screen partition
    graphics::layout(mat = 1,widths = 1,heights = 1)

} else if (which_plot == "peaks_shared") {
    ## Plot peak sharing distribution
out <- hist_shared_peaks(x = x,mcall = mcall,...)
out <- data.frame(rt = out[["time"]], prop = out[["prop"]])
out
}
}
