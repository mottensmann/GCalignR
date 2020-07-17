#' Plot diagnostics for an GCalign Object
#'
#' @description
#' Visualises the aligned data based on four diagnostic plots. One plot shows the distribution of peak numbers per sample in the raw data and after alignment. A second plot gives the distribution of linear shifts that were applied in order to conduct a full alignment of samples with respect to reference. A third sample gives a distribution of the variation in retention times of homologous peaks. The fourth plot shows a frequency distribution of peaks shared among samples.
#'
#' @examples
#' ## GCalign object
#' data("aligned_peak_data")
#'
#' ## All plots are shown by default
#' plot(aligned_peak_data)
#'
#' ## Distribution of peak numbers
#' plot(aligned_peak_data, which_plot = "peak_numbers")
#'
#' ## variation of retention times
#' plot(aligned_peak_data, which_plot = "variation")
#'
#' @param x
#' Object of class GCalign, created with \code{\link{align_chromatograms}}
#'
#' @param which_plot
#' A character defining which plot is created. Options are "shifts", "variation", "peak_numbers" and "peaks_shared". By default all four are created.
#'
#' @param ...
#' Optional arguments passed on to methods. See
#' \code{\link[=graphics]{plot}}, \code{\link[=graphics]{hist}} and \code{\link[=graphics]{barplot}}. Note that optional arguments are not passed on when plotting all figures.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @return
#' Depending on the selected plot a data frame containing the data source of the respective plot is returned. If all plots are created, no output will be returned.
#'
#' @export
#'
plot.GCalign <- function(x, which_plot = c("all","shifts","variation","peak_numbers","peaks_shared"), ...){

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

        # appropriate steps for hist breaks
        if (xmax <= 0.1) bin_size <- 0.01
        if (xmax > 0.1 & xmax < 0.75) bin_size <- 0.025
        if (xmax >= 0.75 & xmax < 5) bin_size <- 0.1
        if (xmax >= 5) bin_size <- NULL

        ## check for optional arguments in the function call, take defaults, if missing
        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Full chromatogram shifts\n(Linear transformation)"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Shift size"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "No. of samples"))
        if (!"breaks" %in% names(mcall) & !is.null(bin_size)) arg_list <- append(arg_list,list(breaks = seq(from = xmin,to = xmax + bin_size, by = bin_size)))
        if (!"freq" %in% names(mcall) & !"frequency" %in% names(mcall)) arg_list <- append(arg_list,list(freq = TRUE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.25))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.25))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#1b9e77"))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))
        # helper to find a good ylim
        p <- as.vector(as.numeric(summary(as.factor(df))))
        p <- max(p)#/sum(p)*100
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))
        if (!"breaks" %in% names(mcall) & !"xaxt" %in% names(mcall) & !is.null(bin_size)) arg_list <- append(arg_list, list(xaxt = "n"))

        x <- do.call(graphics::hist,args = c(list(x = df),arg_list,list(...)))

        if (!"breaks" %in% names(mcall) & !"xaxt" %in% names(mcall) & !is.null(bin_size)) graphics::axis(side = 1, at = x[["mids"]], labels = seq(from = xmin, to = xmax, by = bin_size))

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
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Variation across samples\n(Peak retention time)"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Range (max - min)"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "No. of substances"))
        if (!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(xmin,xmax + 0.01,by = 0.01)))
        if (!"freq" %in% names(mcall) & !"frequency" %in% names(mcall)) arg_list <- append(arg_list,list(freq = TRUE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.25))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.25))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#d95f02"))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"xaxt" %in% names(mcall)) arg_list <- append(arg_list, list(xaxt = "n"))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))
        # helper to find a good ylim
        p <- as.vector(as.numeric(summary(as.factor(df))))
        p <- max(p) #/sum(p)*100
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))
        x <- do.call(graphics::hist,args = c(list(x = df),arg_list,list(...)))
        graphics::axis(side = 1, at = x[["mids"]], labels = seq(xmin, xmax, 0.01))
        return(df)
    }#end function

    bar_peakdistr <- function(x,mcall,...) {
        rt_var_name <- x[["Logfile"]][["Input"]][["Retention_Time"]]
        ## Peaks of all Samples after the alignment
        alg <- (x[["aligned"]][[rt_var_name]])
        ## prior distribution
        pr <- x[["heatmap_input"]][["input_rts"]]
        ## get rid of mean retention time column
        alg <- alg[,2:ncol(alg)]
        pr <- pr[,2:ncol(pr)]
        # transpose pr
        pr <- as.data.frame(t(pr))

        order <- which(names(pr) %in% names(alg)) # Check!

        peak_df <- matrix(0,ncol = 3,nrow = ncol(pr))
        peak_df[,1] <- names(pr)
        peak_df[,2] <- unlist(lapply(1:ncol(pr), function(y) temp <- length(pr[,y][pr[,y] > 0])))
        peak_df[order,3] <- unlist(lapply(1:ncol(alg), function(y) temp <- length(alg[,y][alg[,y] > 0])))
        peak_df <- data.frame(peak_df)
        names(peak_df) <- c("id","pre-aligned","aligned")
        peak_df[["pre-aligned"]] <- as.numeric(as.character(peak_df[["pre-aligned"]]))
        peak_df[["aligned"]] <- as.numeric(as.character(peak_df[["aligned"]]))


        peaks <- peak_df
        peaks[["pre-aligned"]] <- peaks[["pre-aligned"]] - peaks[["aligned"]]
        lablist <- peaks[["id"]]
        peaks <- t(peaks[,c(3,2)])

        ymax <- max(peak_df[,2:3])
        ymin <- min(peaks)

        arg_list <- list()
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = ""))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = ""))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "No. of peaks"))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.25))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.25))
        if (!"cex.names" %in% names(mcall)) {
            lab_thresh <- c(20,30,40,50,60,Inf)
            lab_size <- c(1.2,1.1,0.95,0.85,0.75,0.75)
            samples_size <- ncol(pr)
            # find the matching size
            temp <- which(lab_thresh > samples_size)
            if (min(temp) == 1) {
                label_size <- lab_size[1]
            } else {
                label_size <- lab_size[min(temp) - 1]
            }
        } else {
            label_size <- mcall[["cex.names"]]
        }
        # if (!"cex.names" %in% names(mcall)) arg_list <- append(arg_list,list(cex.names = label_size))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = c("red","blue")))
        # if (!"srt" %in% names(mcall))  arg_list <- append(arg_list,list(srt = 45))
        # if (!"las" %in% names(mcall))  arg_list <- append(arg_list,list(las = 2))
        if (!"names.arg" %in% names(mcall)) arg_list <- append(arg_list,list(names.arg = rep("",each = ncol(peaks))))
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list,list(ylim = c(-3,ymax + 15)))
        # colnames(peaks) <- lablist

        bars <- do.call(graphics::barplot,args = c(list(height = peaks, plot = F),arg_list,list(...)))
        if (!"xlim" %in% names(mcall)) arg_list <- append(arg_list, list(xlim = c(0, max(bars) + 3)))
        bars <- do.call(graphics::barplot,args = c(list(height = peaks),arg_list,list(...)))

        graphics::legend("topleft", rownames(peaks), fill = c("red","blue"), inset = c(-0.009,0), xjust = 0, cex = 0.75, bty = "n")


        if (!"names.arg" %in% names(mcall)) {
            # if ("cex.names" %in% names(mcall)) label_size <- cex.names
            if ("cex.axis" %in% names(mcall)) label_size <- mcall[["cex.axis"]]
            graphics::axis(side = 1, at = bars, labels = lablist, cex.axis = label_size, las = 2)
            #graphics::text(bars, graphics::par("usr")[1], labels = lablist, srt = 90, pos = 1, xpd = TRUE, cex = label_size, offset = 0.4)
        }
        rownames(peak_df) <- peak_df[["id"]]
        return(peak_df)
    }#bar_peak_dist

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
        if (!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Shared substances"))
        if (!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Frequency of samples (%)"))
        if (!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "No. of substances"))
        if (!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(0,101,by = 1)))
        if (!"freq" %in% names(mcall) & !"frequency" %in% names(mcall)) arg_list <- append(arg_list,list(freq = TRUE))
        if (!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis = 1.25))
        if (!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab = 1.25))
        if (!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "grey10"))
        # helper to find a good ylim
        p <- max(as.vector(summary(as.factor(round(peaks_shared$prop)))))
        if (!"ylim" %in% names(mcall)) arg_list <- append(arg_list, list(ylim = c(0,round(p + 5,-1))))
        if (!"right" %in% names(mcall)) arg_list <- append(arg_list,list(right = FALSE))
        if (!"xaxt" %in% names(mcall)) arg_list <- append(arg_list, list(xaxt = "n"))
        if (!"border" %in% names(mcall)) arg_list <- append(arg_list, list(border = "white"))

        x <- do.call(graphics::hist,args = c(list(x = peaks_shared[["prop"]]),arg_list,list(...)))
        graphics::axis(side = 1, at = x[["mids"]][seq(1,101,2)], labels = seq(0, 100, 2))
        return(peaks_shared)
    }#hist_shared_peaks
    # --------------------------------------------------------------------

    if (which_plot == "shifts") {
        out <- hist_linshift(object = x,mcall = mcall,...)
        out <- data.frame(shifts = out)
        #return(out)
    } else if (which_plot == "variation") {
        out <- hist_peakvar(x = x,mcall = mcall,...)
        out <- data.frame(range = as.vector(out), index_sharing = 1:length(as.vector(out)))
        #return(out)
    } else if (which_plot == "peak_numbers") {
        out <- bar_peakdistr(x = x,mcall = mcall,...)
        out <- data.frame(sample = rownames(out), pre_aligned = out[["pre-aligned"]],aligned = out[["aligned"]],row.names = 1:nrow(out))
        #return(out)

    } else if (which_plot == "peaks_shared") {
        ## Plot peak sharing distribution
        out <- hist_shared_peaks(x = x,mcall = mcall,...)
        out <- data.frame(rt = out[["time"]], prop = out[["prop"]])
        #return(out)
    } else if (which_plot == "all")  {
        graphics::layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
        bar_peakdistr(x = x,mcall = mcall)
        hist_linshift(object = x,mcall = mcall)
        hist_peakvar(x = x,mcall = mcall)
        hist_shared_peaks(x = x, mcall = mcall,...)
        ## back to normal screen partition
        graphics::layout(mat = 1,widths = 1,heights = 1)
    }
}
