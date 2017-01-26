#' Visualises peak alignments
#'
#' @description
#' The goal of aligning chromatography peaks is to get the same substance in the same row, although it
#' might have slightly different retention times across samples. This function makes it possible
#' to evaluate the alignment by illustrating how far the retention time of a given peak within a sample
#' deviates from the mean retention time for this substance across all samples after alignment. In other words,
#' given that the mean retention time is the 'real' retention time of a substance, the heatmap
#' shows how far a given peak deviates from it. Two types of heatmaps are available. A binary heatmap allows to determine
#' if the retention time of single samples deviates by more than a user defined threshold from the mean. Optionally, a discrete heatmap allows to check deviations quantitatively.
#'
#' @param object
#' Object of class "GCalign", the output of a call to \link{align_chromatograms}.
#'
#' @param algorithm_step
#' Character indicating which step of the algorithm is plotted. Either \strong{pre_alignment}, \strong{linear_shifted} or \strong{aligned} specifiying the raw, linearly shifted or aligned data respectively.
#' Default is the heatmap after alignment.
#'
#' @param substance_subset
#' Vector containing indices of substances (ordered in ascending order of retention times) to plot. By default \code{NULL} indicates that all substances are plotted.
#'
#' @param legend_type
#' Character specifying the type of colourbar as \strong{discrete} (i.e retention times are classified as deviating or not) or \strong{gradient} (i.e deviations are presented on a fine scale).
#'
#' @param  samples_subset
#' Vector indicating which samples are plotted on the heatmap.
#' Either a numeric vector of indices (order in the input) or a vector of sample names.
#'
#' @param  type
#' Character specifying whether a \strong{'binary'} heatmap or a heatmap of \strong{'discrete'}
#' deviations is plotted.
#'
#' @param threshold
#' Numerical value denoting the threshold above which the deviation of individual peak retention times
#' from the mean retention time of the respective substance is highlighted in \emph{binary} heatmaps.
#'
#'@param label_size
#' Integer determining the size of labels on y and x axis. By default the label_size is calculated (beta!) to compromise between readibility and messines due to a potentially large number of substances and samples. Note: Labels for substances on the x axis are only possible if a maximum of 150 substances are plotted.
#'
#'@param show_legend
#' Logical determining whether a legend is included or not. Default is TRUE.
#'
#'
#'@param main_title
#' Character argument used a title for the plot. If not specified, titles are generated automatically.
#'
#' @param label
#' Logical determining whether labels for samples and substances are used. By default sample names are written on the vertical axis, but in large data sets its less readible and may be used effectively only for subsets of the data.
#'
#'@return
#' object of class "ggplot"
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import ggplot2 RColorBrewer grDevices
#'
#' @examples
#'
#'  ## Default settings: The final output is plotted
#'  gc_heatmap(aligned_peak_data, algorithm_step="aligned")
#'
#'  ## Plot the input data
#'  gc_heatmap(aligned_peak_data,algorithm_step="pre_alignment")
#'
#'  ## Plot a subset of the first 50 scored substances
#'  gc_heatmap(aligned_peak_data,algorithm_step="aligned",substance_subset = 1:50)
#'
#'  ## Plot specific samples, apply a stricter threshold
#'  gc_heatmap(aligned_peak_data,samples_subset = c("M2","P7","M13","P13"),threshold=0.02)
#'
#' @export
#'
gc_heatmap <- function(object, algorithm_step = c('aligned','linear_shifted','pre_alignment'),
    substance_subset = NULL, legend_type = c('legend','colourbar'), samples_subset = NULL,
    type = c("binary","discrete"), threshold = 0.05, label_size = NULL, show_legend = TRUE,
    main_title = NULL, label = TRUE) {

    algorithm_step <- match.arg(algorithm_step)
    if (algorithm_step == "aligned") algorithm_step <- "aligned_rts"
    if (algorithm_step == "linear_shifted") algorithm_step <- "linear_transformed_rts"
    if (algorithm_step == "pre_alignment") algorithm_step <- "input_rts"

    type <- match.arg(type)

    legend_type <- match.arg(legend_type)



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
heat_matrix['substance'] <- as.factor(round(as.numeric(as.character(heat_matrix[['substance']])),digits = 2))

# Plot the heatmap
    if (type == "binary") {
        # Case that no deviations are present, i.e. perfect alignment
        if (max(heat_matrix['diff'],na.rm = T) == 0) {
            hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id',fill = 'diff'),colour = "Blue")
            hm <- hm + geom_tile(color = "transparent", size = 0.001)
            hm <- hm + scale_fill_gradientn(colours = 'blue',na.value = "white")
            hm <- hm + labs(x = "Substances", y = "Samples", title = ifelse(is.null(main_title),paste("No deviations exceeding a threshold of",as.character(threshold)),main_title))
            hm <- hm + guides(fill = FALSE)
        # Usual Case, at leat some samples deviate at certain retention times
        } else {
            hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id',fill = 'diff'))
            hm <- hm + geom_tile(color = "transparent", size = 0.001)
            hm <- hm + scale_fill_continuous(low = "#a6cee3",high = "#b2182b",breaks = c(0,1),na.value = "white", guide = 'legend',name = paste('Deviation\n','>',as.character(threshold)),labels = c('NO','YES'))
            hm <- hm + labs(x = "Substances", y = "Samples", title = ifelse(is.null(main_title),paste("Deviation from substance mean retention time\n(Threshold = ",as.character(threshold),")"),main_title))
        }
        # type == continuos
    } else {
        # Spectral colour scheme
        myPalette <-  colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
        r <- c(min(heat_matrix[["diff"]],na.rm = T),max(heat_matrix[["diff"]],na.rm = T))

        # Take 10 colours form the spectral scheme
        # colourset <- myPalette(10)
        # remove middle ones for higher contrast
        # colourset <- colourset[c(1:4,7:10)]

        hm <- ggplot(heat_matrix, aes_string(x = 'substance', y = 'id', fill = 'diff'))
        hm <- hm + ggplot2::geom_tile(color = "white", size = 0.01)
        hm <- hm + scale_fill_gradientn(colours = myPalette(10),guide = "legend",name = 'Deviation',na.value = "white", limits = c(-round(max(abs(r)),2) - 0.01,round(max(abs(r)),2) + 0.01)
        )
        hm <- hm + labs(x = "Substances", y = "Samples", title = ifelse(is.null(main_title),"Variation of retention times",main_title))
    }
    hm <- hm + theme(plot.title = element_text(hjust = 0.5,vjust = 1,size = 10,face = 'bold'))
    hm <- hm + theme(axis.title.x = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.x = element_blank(),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
                     axis.text.y = element_text(size = label_size,hjust = 0.5))
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
if ((!is.null(substance_subset) & ncol(rt_df) < 151) || ncol(rt_df) < 151 ) {
    hm <- hm + theme(axis.text.x = element_text(size = label_size, hjust = 0.5,angle = 90),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                     axis.text.y = element_text(size = label_size,hjust = 0.5))
}
if (!show_legend) {hm <- hm + theme(legend.position = "none")}
if (!label) {hm <- hm + theme(axis.text = element_blank())}
# return the ggplot object
return(hm)
}
