#' Plot Diagonstics for an Gcalign Object
#'
#' @description
#' Three plots are currently available: One plot visualises the distribution of linear shifts
#' that were applied to align chromatograms to a reference before aligning individual peaks.
#' A second plot illustrates the remaining variation of retention times on the level of individual
#' peaks by plotting the distribution of retention time ranges. The third plots shows a distribution
#' of peak numbers after aligning the chromatograms.
#'
#' @examples
#' ## All three plots
#' plot(aligned_peak_data)
#'
#' ## Distribution of Peaks
#' plot(aligned_peak_data,which_plot="Peak_Counts")
#'
#' @param x \code{GCalign} object, created by \code{\link{align_chromatograms}}
#'
#' @param which_plot
#' character string indicating which plot is returned. Available are
#' \strong{"Linear_Shifts"} a histogram of linear adjustments undertaken in aligning chromatograms,
#' \strong{"Peak_Range"} a histogram summarising the range of retention times for every peak defined
#' by the difference between minimum and maximum retention times respectively. The third option
#' is \strong{"Peak_Counts"} plotting a barchart of the number of peaks per sample. By default all
#' three plots are returned as subplots of one figure.
#'
#' @param ...
#' optional arguments passed on to methods. See
#' \code{\link[graphics]{plot}}, \code{\link[graphics]{hist}} and \code{\link[graphics]{barplot}}.
#' Please Note that optional arguments are currently not passed on when plotting all figures.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @export
#'
plot.GCalign <- function(x,which_plot=c("All","Linear_Shifts","Peak_Range","Peak_Counts","Peak_Sharing"),...){

    # initialising by picking arguments from the function call
    mcall = as.list(match.call())[-1L]
    which_plot <- match.arg(which_plot) # allow partial matching


    # Define internal functions
# -------------------------------------------------------------------

    hist_linshift <- function(object,mcall,...){

         xl <- object[["Logfile"]][["LinearShift"]]["shift"] # Get the applied shifts

        df <- as.vector(unlist(as.vector(xl))) # steps of linear shifts and their frequency
        xmax <- object[["Logfile"]][["Call"]][["max_linear_shift"]]
        xmin <- -xmax

        # check for optional arguments in the function call, take defaults, if missing
        arg_list <- list()
        if(!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Linear Transformation"))
        if(!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Shift"))
        if(!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Frequency [%]"))
        if(!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks = seq(from=xmin,to=xmax,by=0.01)))
        if(!("freq") %in% names(mcall)||!("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = FALSE))
        if(!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis=1.5))
        if(!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab=1.5))
        if(!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#ffffbf"))


        do.call(graphics::hist,args=c(list(x=df),arg_list,...))
    }#end hist_linshift

    hist_peakvar <- function(x,mcall,...){
    MinMax <- function(x){ # Estimate the range of retention times per substance, they should be no overlapp
        temp <- matrix(NA,1,2)
        colnames(temp) <- c("range","row")
        data <- temp[0,]
        for(i in 1:ncol(x)){
            data<-rbind(data,cbind(abs(diff(range(x[,i][x[,i]>0],na.rm = TRUE))),i))
        }
        return(as.data.frame(data))
    }

        df <- MinMax(x[["heatmap_input"]][["aligned_rts"]][,-1]) # Range of RTs aligned
        df <- unlist(df["range"]) # Formatting

        xmax <- max(df)
        xmin <- min(df)

        arg_list <- list()
        if(!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Retention Time Range\nper Peak"))
        if(!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Range [max - min]"))
        if(!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Frequency [%]"))
        if(!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks =seq(xmin,xmax,by=0.01)))
        if(!("freq") %in% names(mcall)||!("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = FALSE))
        if(!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis=1.5))
        if(!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab=1.5))
        if(!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#91bfdb"))

        do.call(graphics::hist,args=c(list(x=df),arg_list,...))

        # xfit <- seq(xmin,xmax,by=0.01)
        # yfit <- dnorm(xfit,mean=mean(df),sd=sd(df))
        # #yfit <- yfit*diff(hist[["mids"]][1:2]*length(df))
        # lines(xfit,yfit,col="black",lwd=2)


    }

    bar_peakdistr <- function(x,mcall,...){
        rt_var_name <- x[["Logfile"]][["Input"]][["Retention_Time"]]
        conc_var_name <- x[["Logfile"]][["Input"]][["Concentration"]]

        data <- (x[["aligned"]][[rt_var_name]]) # Peaks of All Samples
        data <- data[,2:ncol(data)] # get rid of mean retention time column

        peak_df <- matrix(NA,ncol = 2,nrow = length(data))
        peak_df[,1] <- names(data)
        peak_df[,2] <- unlist(lapply(1:ncol(data), function(y) temp <- length(data[,y][data[,y]>0])))
        peak_df <- data.frame(peak_df)
        names(peak_df) <- c("ID","Peaks")
        peak_df$Peaks <- as.numeric(as.character(peak_df$Peaks))

        peaks <- peak_df[["Peaks"]]
        names(peaks) <- peak_df[["ID"]]

        ymax <- max(peaks)
        ymin <- min(peaks)

        arg_list <- list()
        if(!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Peak Counts\n after running GCalignR"))
        if(!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = ""))
        if(!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Number of Peaks"))
        if(!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis=1.5))
        if(!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab=1.5))
        if(!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#fc8d59"))
        if(!"srt"%in% names(mcall))  arg_list <- append(arg_list,list(srt = 45))
        if(!"las"%in% names(mcall))  arg_list <- append(arg_list,list(las = 2))
        if(!"names.arg" %in% names(mcall)) arg_list <- append(arg_list,list(names.arg=names(peaks)))
        if(!"ylim" %in% names(mcall)) arg_list <- append(arg_list,list(ylim=c(0,ymax+5)))

        bars <- do.call(graphics::barplot,args=c(list(height=peaks),arg_list,...))
        graphics::text(x=bars,y=peaks+2,labels = as.character(peaks),cex = 0.9)

    }
    hist_shared_peaks <- function(x,mcall,...){
        df <- x[["heatmap_input"]][["aligned_rts"]][-1] # check that this is up to date!
        shared_substances <- data.frame(time=round(as.numeric(names(df)),3),prop=unlist(lapply(1:ncol(df), function(col){
            sub <- df[,col] # all entries of a column
            N <- length(sub[sub>0])
            N/length(sub)*100
        })))


        # xmax <- max(shared_substances["prop"])
        # xmin <- min(shared_substances["prop"])

        arg_list <- list()
        if(!"main" %in% names(mcall)) arg_list <- append(arg_list,list(main = "Shared Peaks among Samples"))
        if(!"xlab" %in% names(mcall)) arg_list <- append(arg_list,list(xlab = "Number of Substances"))
        if(!"ylab" %in% names(mcall)) arg_list <- append(arg_list,list(ylab = "Proportion of Samples [%]"))
        if(!"breaks" %in% names(mcall)) arg_list <- append(arg_list,list(breaks =seq(0,100,by=5)))
        if(!("freq") %in% names(mcall)||!("frequency") %in% names(mcall)) arg_list <- append(arg_list,list(freq = FALSE))
        if(!"cex.axis" %in% names(mcall)) arg_list <- append(arg_list,list(cex.axis=1.5))
        if(!"cex.lab" %in% names(mcall)) arg_list <- append(arg_list,list(cex.lab=1.5))
        if(!"col" %in% names(mcall))  arg_list <- append(arg_list,list(col = "#91bfdb"))

        do.call(graphics::hist,args=c(list(x=shared_substances[["prop"]]),arg_list,...))

    }
# --------------------------------------------------------------------

if(which_plot=="Linear_Shifts"){
hist_linshift(object = x,mcall = mcall,...)

}else if (which_plot=="Peak_Range"){
    hist_peakvar(x = x,mcall = mcall,...)

} else if(which_plot=="Peak_Counts"){
    bar_peakdistr(x = x,mcall = mcall,...)

} else if(which_plot=="All"){
    graphics::layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    bar_peakdistr(x = x,mcall = mcall)
    hist_linshift(object = x,mcall = mcall)
    hist_peakvar(x = x,mcall = mcall)
    graphics::layout(mat = 1,widths = 1,heights = 1)#back to normal screen partition

} else{
    # Plot peak sharing distribution
hist_shared_peaks(x=x,mcall=mcall,...)
}
}


