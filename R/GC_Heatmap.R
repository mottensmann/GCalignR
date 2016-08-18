#' gc_heatmap visualises the goodness of a chromatogram alignments
#'
#' @description Visualises the deviation of single samples from the whole population by comparing
#' the actual retention time against the mean of all other samples containing the same
#' substance. Two types of heatmaps are available. A binary heatmap allows to determine
#' if single samples within each chromatogram are correctly assigned to a certain substance,
#' by setting a fixed threshold of allowed deviations from the mean. The optional discrete
#' heatmap allows to check the deviations quantitatively.
#'
#'
#' @param object object of class "GCalign", the output of a call to \link{align_chromatograms}.
#'
#' @param algorithm_step
#' character indicating which step of the algorithm is plotted. Either \strong{input_rts}, \strong{linear_transformed_rts} or \strong{aligned_rts}.
#'
#' @param substance_subset
#' Vector containing indices of substances (ordered in ascending order of retention times) to plot. By default \code{NULL} indicating all substances are plotted.
#'
#' @param guide
#' Character string, selects the type of colourbar as discrete (i.e \strong{'legend'})
#' or gradient (i.e \strong{'colourbar})
#'
#' @param  samples_subset
#' Vector indicating which samples are plotted on the heatmap.
#' Either a numeric vector of indices (order in the input) or a vector of sample names
#'
#' @param  type
#' Character specifying whether a binary heatmap or a heatmap of continous
#' deviations is plotted.
#'
#' @param threshold
#' Decimal indicating the threshold deviation of individual peak retention times
#' from the mean retention time of the respective peak across all samples.
#'
#'@param label_size
#'Determines the size of labels on y and x axis. By default the label_size is calculated (beta!) to compromise between readibility and messines due to a potentially large number of substances and samples. Note: Label for substances on the x axis are only plotted if a subset of substances was selected, or less than xxx substances are present in the data.
#'
#' @return
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
#'  gc_heatmap(aligned_peak_data, algorithm_step="aligned_rts")
#'
#'  ## Plot the input data
#'  gc_heatmap(aligned_peak_data,algorithm_step="input_rts")
#'
#'  ## Plot a subset of the first 50 scored substances
#'  gc_heatmap(aligned_peak_data,algorithm_step="aligned_rts",substance_subset = 1:50)
#'
#'  ## Plot specific samples, apply a stricter threshold
#'  gc_heatmap(aligned_peak_data,samples_subset = c("M2","P7","M13","P13"),threshold=0.02)
#'
#' @export
#'

gc_heatmap <-function(object,algorithm_step=c('aligned_rts','linear_transformed_rts','input_rts'),substance_subset=NULL,guide=c('legend','colourbar'),samples_subset=NULL,type=c("binary","continous"),threshold=0.05,label_size=NULL){

## Grip called parameters
    algorithm_step <- match.arg(algorithm_step)
    type <- match.arg(type)
    guide <- match.arg(guide)

## Get the retention time matrix for the selected step
    rt_df <- object[['heatmap_input']][[algorithm_step]]
    rt_df[,'id'] <- as.character(rt_df[,'id'])

## Try to estimate a suitable label_size
if(is.null(label_size)){
    lab_thresh <- c(20,40,60,80,100,120,140,999) # should be coded saver!
    lab_size <- c(12,10,8,8,6,5.5,4,4)
    samples_size <- nrow(rt_df)
    ## find the matching size
    temp <- which(lab_thresh > samples_size)
    label_size <- lab_size[min(temp)-1]
}
## Get the initial order, already coded in align_chromatograms
    # rt_df <- rt_df[match(as.character(object[["heatmap_input"]][["input_rts"]][,1]),as.character(rt_df[,1])),]

## Select subsets if specified
    if(!is.null(substance_subset)){ # substances
        rt_df <- rt_df[,c(1,substance_subset+1)] # always keep id !, therefore +1
    }

    if(!is.null(samples_subset)){ # samples
        if(is.character(samples_subset)){
            rt_df <- rt_df[rt_df[,1]%in%samples_subset,]
        }else if(is.numeric(samples_subset)){
            rt_df <- rt_df[samples_subset,]
        }
    }
## Create a data frame for plotting with ggolot

    heat_matrix <- reshape2::melt(data = rt_df,id.vars ='id')
    names(heat_matrix) <- c('id','substance','rt')
    heat_matrix[,'substance'] <- as.numeric(as.character(heat_matrix[,'substance']))

## Calculate the deviation of each peak from its substance mean
    heat_matrix[,'diff'] <- (as.numeric(heat_matrix[,'rt']) - heat_matrix[,'substance'])
    heat_matrix['diff'][heat_matrix['rt']==0] <- 0 # zero mean no substance is present

    heat_matrix[,'id'] <- ordered( heat_matrix[,'id'], levels = as.factor(rt_df[,'id']))
    heat_matrix[,'substance'] <- ordered( heat_matrix[,'substance'], levels = as.factor(colnames(rt_df)[2:ncol(rt_df)]))



## If binary heatmap was selected, code violoation at the level of the threshold by 0/1
    if(type=="binary"){

        heat_matrix['diff'][abs(heat_matrix['diff'])>threshold] <- 1 # Deviates
        heat_matrix['diff'][abs(heat_matrix['diff'])<threshold] <- 0 # Is okay
        heat_matrix['diff'][heat_matrix['rt']==0] <- NA
    }
## Code absence by NA
heat_matrix['diff'][heat_matrix['rt']==0] <- NA

## Simplify substance names
heat_matrix['substance'] <- as.factor(round(as.numeric(as.character(heat_matrix[['substance']])),digits = 2))

## Plot the heatmap

    if(type=="binary"){
        if(max(heat_matrix['diff'],na.rm = T)==0){ # Zeros indicates no deviations of any substance
            hm <- ggplot(heat_matrix, aes_string(x='substance', y='id',fill='diff'),colour="Blue")
            hm <- hm + geom_tile(color="transparent", size=0.001)
            hm <- hm + scale_fill_gradientn(colours = 'blue',na.value = "white")
            hm <- hm + labs(x=NULL, y=NULL, title=paste("No deviations exceeding a threshold of",as.character(threshold)))
            hm <- hm + guides(fill=FALSE)


## Case is not possible, therefore exluded
        # }else if(min(heat_matrix['diff'],na.rm = T)==1){ # Really bad alignment
        #     hm <- ggplot(heat_matrix, aes_string(x='substance', y='id',fill=diff),colour="red")
        #     hm <- hm + geom_tile(color="transparent", size=0.001)
        #     hm <- hm + scale_fill_gradientn(colours = 'red',na.value = "white")
        #     hm <- hm + labs(x=NULL, y=NULL, title=paste("Alignment failed, all substances show deviation from threshold of",as.character(threshold)))
        #     hm <- hm + guides(fill=FALSE)

        }else{ # At leat some samples deviate at certain retention times
            hm <- ggplot(heat_matrix, aes_string(x='substance', y='id',fill='diff'))
            hm <- hm + geom_tile(color="transparent", size=0.001)
            hm <- hm + scale_fill_continuous(low = "#a6cee3",high = "#b2182b",breaks=c(0,1),na.value = "white",
                                             guide = 'legend',name=paste('Deviation\n','>',as.character(threshold)),labels=c('NO','YES'))
            hm <- hm + labs(x=NULL, y=NULL, title=paste("Deviations of retention times at a threshold of",as.character(threshold)))
        }
    }else{
        myPalette <-  colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral"))) # Spectral colours
        colourset <- myPalette(10) # Take 10 colours form the spectral scheme
        colourset <- colourset[c(1:4,7:10)] # remove middle ones for higher contrast to background
        hm <- ggplot(heat_matrix, aes_string(x='substance', y='id', fill='diff'))
        hm <- hm + geom_tile(color="white", size=0.01)
        hm <- hm + scale_fill_gradientn(colours = myPalette(10),guide = guide,name='Deviation',na.value = "white")
        hm <- hm + labs(x=NULL, y=NULL, title="Deviation of individual peak retentention times from the substance mean")

    }
    hm <- hm + theme(plot.title=element_text(hjust=0.5,vjust=1,size = 16,face = 'bold'))
    hm <- hm + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
                     axis.text.y=element_text(size = label_size,hjust = 0.5))
    #hm <- hm + coord_equal(ncol(rt_df)/nrow(rt_df))

hm <- hm + theme(plot.background = element_rect(fill = "grey95"))

y <- 1:nrow(rt_df)+.5
x <- rep(0,nrow(rt_df))
yend <- 1:nrow(rt_df)+.5
xend <- rep(ncol(rt_df),nrow(rt_df))
my.lines<-data.frame(y=y,x=x,xend=xend, yend=yend)

hm <- hm + geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend),color="grey", size=0.35,show.legend = FALSE, inherit.aes=F)

y <- rep(0,ncol(rt_df))
x <- 1:ncol(rt_df)+.5
xend <- 1:ncol(rt_df)+.5
yend <- rep(nrow(rt_df)+.5,ncol(rt_df))
my.lines <- data.frame(y=y,x=x, xend=xend,yend=yend)

hm <- hm + geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend),color="grey", size=0.35,show.legend = FALSE, inherit.aes=F)



## If subsets are selected, allow to plot the substance labels for better idetenfication
if(!is.null(substance_subset)){
    hm <- hm + theme(axis.title.x=element_blank(),
                     axis.text.x=element_text(size = label_size,hjust = 0.5,angle = 90),
                     axis.ticks.y = element_line(size = 0.3, colour = "grey60"),
                     axis.ticks.x = element_line(size = 0.3, colour = "grey60"),
                     axis.text.y=element_text(size = label_size,hjust = 0.5))
}

return(hm)

}


