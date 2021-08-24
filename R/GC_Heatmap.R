#' Visualises peak alignments in form of a heatmap
#'
#' @description
#' The goal of aligning peaks is to match homologous peaks that are thought to represent homologous substances in the same row across samples, although peaks have slightly different retention times across samples. This function makes it possible to evaluate the alignment quickly by inspecting the (i) distribution of peaks across samples, (ii) the variation for each homologous peak (column) as well as (iii) patterns that might hint at splitting peaks across rows. The mean retention time per homologous peak is here defined as the "true"  retention time and deviations of individual peaks can be seen by a large deviation in the retention time to the mean. Subsetting of the retention time range (i.e. selecting peaks by the mean retention time) and samples (by name or by position) allow to quickly inspect regions of interest. Two types of heatmaps are available, a binary heatmap allows to determine if the retention time of single samples deviates by more than the user defined threshold from the mean. Optionally, a discrete heatmap allows to check deviations quantitatively. Large deviation can have multiple reasons. The most likely explanation is given by the fact that adjacent rows were merged as specified by the value \code{min_diff_peak2peak} in \code{\link{align_chromatograms}}. Here clear cases, in which peaks of multiple samples have been grouped in either of two or more rows can be merged and cause relatively high variation in peak retention times.
#'
#' @param object
#' Object of class "GCalign", the output of a call to \code{\link{align_chromatograms}}.
#'
#' @param algorithm_step
#' Character indicating which step of the algorithm is plotted. Either "input", "shifted" or "aligned" specifying the raw, linearly shifted or aligned data respectively. Default is the heatmap for the aligned dataset.
#'
#' @param substance_subset
#' A vector of integers containing indices of substances in ascending order of retention times to plot.
#'
#' @param legend_type
#' A character specifying how to present deviations of retention times from the mean. Either in form of discrete steps or on a gradient scale using 'legend' or 'colourbar' respectively. Changes are only possible when \code{type = "discrete"}
#'
#' @param  samples_subset
#' A vector indicating which samples are plotted on the heatmap by giving either indices or names of samples.
#'
#' @param  type
#' A character specifying whether a deviations of retention times are shown 'binary' (i.e. in comparison to the threshold value) or on a 'discrete' scale with respect to the mean retention time.
#'
#' @param threshold
#' A numeric value denoting the threshold above which the deviation of individual peak retention times from the mean retention time of the respective substance are highlighted in heatmaps. By default, the value of parameter \code{max_diff_peak2mean} (see \code{\link{align_chromatograms}}) that was used in aligning the data is used.
#'
#'@param label_size
#' An integer determining the size of labels on y and x axis. By default a fitting label_size is calculate (beta!) to compromise between readability and messiness due to a potentially large number of substances and samples.
#'
#'@param show_legend
#' Boolean determining whether a legend is included or not.
#'
#'@param main_title
#' Character giving the title of the heatmap. If not specified, titles are generated automatically.
#'
#' @param label
#' Character determining if labels are shown on axes. Depending on the number of peaks and/or samples, labels are difficult to read. Use subsets instead. Possible option are "xy", "x", "y" or "none"
#'
#'@return
#' object of class "ggplot"
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import ggplot2
#'
#' @examples
#'
#'  ## aligned gc-dataset
#'  data("aligned_peak_data")
#'  ## Default settings: The final output is plotted
#'  gc_heatmap(aligned_peak_data, algorithm_step = "aligned")
#'
#'  ## Plot the input data
#'  gc_heatmap(aligned_peak_data,algorithm_step = "input")
#'
#'  ## Plot a subset of the first 50 scored substances
#'  gc_heatmap(aligned_peak_data,algorithm_step="aligned",substance_subset = 1:50)
#'
#'  ## Plot specific samples, apply a stricter threshold
#'  gc_heatmap(aligned_peak_data,samples_subset = c("M2","P7","M13","P13"),threshold = 0.02)
#'
#' @export
#'
gc_heatmap <- function(object = NULL,
                       algorithm_step = c('aligned','shifted','input'),
                       substance_subset = NULL, legend_type = c('legend','colourbar'),
                       samples_subset = NULL,
                       type = c("binary","discrete"),
                       threshold = NULL,
                       label_size = NULL,
                       show_legend = TRUE,
                       main_title = NULL,
                       label = c("y","xy","x","none")) {

    # removed from @import: RColorBrewer grDevices
    algorithm_step <- match.arg(algorithm_step)
    if (algorithm_step == "aligned") algorithm_step <- "aligned_rts"
    if (algorithm_step == "shifted") algorithm_step <- "linear_transformed_rts"
    if (algorithm_step == "input") algorithm_step <- "input_rts"

    if (is.null(threshold)) threshold <- object[["Logfile"]][["Call"]][["max_diff_peak2mean"]]
    type <- match.arg(type)
    legend_type <- match.arg(legend_type)
    label <- match.arg(label)

    # Get the retention time matrix for the selected step
    rt_df <- object[['heatmap_input']][[algorithm_step]]
    rt_df[,'id'] <- as.character(rt_df[,'id'])

    # Try to estimate a suitable label_size
    if (is.null(label_size)) {
        lab_thresh <- c(20,40,60,80,100,120,140, Inf)
        lab_size <- c(12,10,8,8,6,5,4,4)
        samples_size <- nrow(rt_df)
        # find the matching size
        temp <- which(lab_thresh > samples_size)
        if (min(temp) == 1) {
            label_size <- lab_size[1]
        } else {
            label_size <- lab_size[min(temp) - 1]
        }
    }

    # Substance subsetting
    if (!is.null(substance_subset)) {
        # always keep id column, therefore + 1
        rt_df <- rt_df[,c(1,substance_subset + 1)]
    }
    # samples subsetting
    if (!is.null(samples_subset)) {
        if (is.character(samples_subset)) {
            rt_df <- rt_df[rt_df[,1] %in% samples_subset,]
        } else if (is.numeric(samples_subset)) {
            rt_df <- rt_df[samples_subset,]
        }
    }
    # Create a data frame for plotting with ggolot
    heat_matrix <- reshape2::melt(data = rt_df,id.vars = 'id')
    names(heat_matrix) <- c('id','substance','rt')
    heat_matrix[,'substance'] <- as.numeric(as.character(heat_matrix[,'substance']))

    # Calculate the deviation of each peak from its substance mean
    heat_matrix[,'diff'] <- (as.numeric(heat_matrix[,'rt']) - heat_matrix[,'substance'])
    # zero means no substance is present
    heat_matrix['diff'][heat_matrix['rt'] == 0] <- 0
    heat_matrix[,'id'] <- ordered( heat_matrix[,'id'], levels = as.factor(rt_df[,'id']))
    heat_matrix[,'substance'] <- ordered( heat_matrix[,'substance'], levels = as.factor(colnames(rt_df)[2:ncol(rt_df)]))

    # If binary heatmap was selected, code violoation at the level of the threshold by 0/1
    if (type == "binary") {
        # Deviation coded by 1
        heat_matrix['diff'][abs(heat_matrix['diff']) > threshold] <- 1
        # No deviation coded by 0
        heat_matrix['diff'][abs(heat_matrix['diff']) < threshold] <- 0
    }
    # Code absence by NA
    heat_matrix['diff'][heat_matrix['rt'] == 0] <- NA

    # Simplify substance names
    heat_matrix['substance'] <- as.factor(as.numeric(as.character(heat_matrix[['substance']])))
    #heat_matrix['substance'] <- as.factor(round(as.numeric(as.character(heat_matrix[['substance']])), digits = 2))

    # Plot the heatmap
    if (type == "binary") {
        # Case that no deviations are present, i.e. perfect alignment
        if (max(heat_matrix['diff'],na.rm = T) == 0) {
            hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id',fill = 'diff'),colour = "Blue")
            hm <- hm + geom_tile(color = "transparent", size = 0.001)
            hm <- hm + scale_fill_gradientn(colours = 'blue',na.value = "white")
            hm <- hm + labs(x = "substance", y = "sample", title = ifelse(is.null(main_title),paste("No deviations exceeding a threshold of",as.character(threshold)),main_title))
            hm <- hm + guides(fill = FALSE)
            # Usual Case, at leat some samples deviate at certain retention times
        } else if (all(heat_matrix['diff'][!is.na(heat_matrix['diff'])] == 1)) {
            hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id',fill = 'diff'),colour = "Blue")
            hm <- hm + geom_tile(color = "transparent", size = 0.001)
            hm <- hm + scale_fill_gradientn(colours = 'red',na.value = "white")
            hm <- hm + labs(x = "substance", y = "sample", title = ifelse(is.null(main_title),paste("Deviations exceeding threshold of",as.character(threshold)),main_title))
            hm <- hm + guides(fill = FALSE)
        } else {
        hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id',fill = 'diff'))
        hm <- hm + geom_tile(color = "transparent", size = 0.001)
        hm <- hm + scale_fill_continuous(low = "blue",high = "red",breaks = c(0,1),na.value = "white", guide = 'legend',name = paste('Deviation\n','>',as.character(threshold)),labels = c('NO','YES'))
        hm <- hm + labs(x = "substance", y = "sample", title = ifelse(is.null(main_title),paste("Deviation from substance mean retention time\n(Threshold = ",as.character(threshold),")"),main_title))
    }
    # type == continuous
} else {
    col_pal <- c("#5E4FA2","#378EBA","#75C8A4","#BEE4A0","#F1F9A9","#FEEDA2", "#FDBE6F","#F67B49","#D8434D","#9E0142")
    r <- c(min(heat_matrix[["diff"]],na.rm = T),max(heat_matrix[["diff"]],na.rm = T))


    hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id', fill = 'diff'))
    hm <- hm + geom_tile(color = "white", size = 0.01)
    hm <- hm + scale_fill_gradientn(colours = col_pal,guide = "legend",name = 'Deviation',na.value = "white", limits = c(-round(max(abs(r)),2) - 0.01,round(max(abs(r)),2) + 0.01)
    )
    hm <- hm + labs(x = "substance", y = "sample", title = ifelse(is.null(main_title),"Variation of retention times",main_title))
}
hm <- hm + theme(plot.title = element_text(hjust = 0.5,vjust = 1,size = 10,face = 'bold'))


theme(axis.text.x = element_text(size = label_size, hjust = 0.5,angle = 90),
      axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
      axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
      axis.text.y = element_text(size = label_size,hjust = 0.5))



if (label == "xy") {
    hm <- hm + theme(axis.title.x = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.x = element_text(size = label_size, vjust = 0.5,angle = 90),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                     axis.text.y = element_text(size = label_size,hjust = 0.5))

} else if (label == "y") {
    hm <- hm + theme(axis.title.x = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.x = element_blank(),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                     axis.text.y = element_text(size = label_size,hjust = 0.5))

} else if (label == "x") {
    hm <- hm + theme(axis.title.x = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.x = element_text(size = label_size, vjust = 0.5,angle = 90),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                     axis.text.y = element_blank())

} else if (label == "none") {
    hm <- hm + theme(axis.title.x = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.x = element_blank(),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
                     axis.text.y = element_blank())

}
hm <- hm + theme(plot.background = element_rect(fill = "grey95"))

# Scoping issues reuire the following loop way to define the data frame
y <- 1:nrow(rt_df) + 0.5
x <- rep(0,nrow(rt_df))
yend <- 1:nrow(rt_df) + 0.5
xend <- rep(ncol(rt_df),nrow(rt_df))
my.lines <- data.frame(y = y,x = x,xend = xend, yend = yend)

hm <- hm + geom_segment(data = my.lines, aes(x,y,xend = xend, yend = yend),color = "grey", size = 0.35,show.legend = FALSE, inherit.aes = F)

y <- rep(0,ncol(rt_df))
x <- 1:ncol(rt_df) + 0.5
xend <- 1:ncol(rt_df) + 0.5
yend <- rep(nrow(rt_df) + 0.5,ncol(rt_df))
my.lines <- data.frame(y = y,x = x, xend = xend,yend = yend)

hm <- hm + geom_segment(data = my.lines, aes(x,y,xend = xend, yend = yend),color = "grey", size = 0.35,show.legend = FALSE, inherit.aes = F)

# If subsets are selected, allow to plot the substance labels for better idetenfication
if (is.null(label)) {
    if ((!is.null(substance_subset) & ncol(rt_df) < 151) || ncol(rt_df) < 151 ) {
        hm <- hm + theme(axis.text.x = element_text(size = label_size, vjust = 0.5,angle = 90),
                         axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                         axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                         axis.text.y = element_text(size = label_size,hjust = 0.5))
    }
}

if (!show_legend) hm <- hm + theme(legend.position = "none")
# return the ggplot object
return(hm)
}
