#' GC_Heatmap visualises the goodness of the chromatogram alignment
#'
#' @description Visualises the deviation of single samples from the whole population by comparing
#'          the actual retention time against the mean of all other samples containing the same
#'          substance. Two types of heatmaps are available. A binary heatmap allows to determine
#'          if single samples within each chromatogram are correctly assigned to a certain substance,
#'          by setting a fixed threshold of allowed deviations from the mean. The optional discrete
#'          heatmap allows to check the deviations quantitatively.
#'
#' @param GcOut \code{data.frame} representing a matrix of retention times,
#'          where samples are ordered in rows
#'
#' @param algorithm_step \code{character} indicating which step of the algorithm is plotted. Either
#'          \code{rt_raw}, \code{rt_linear} or \code{rt_aligned}
#'
#' @param substance_subset \code{vector} containing indices of substances (i.e. rows) to plot
#'          By default \code{NULL} indicating all substances are plotted
#'
#' @param guide \code{character} indicating type of colourbar as discrete (i.e 'legend')
#'          or gradient (i.e 'colourbar)
#'
#' @param  samples_subset \code{vector} indicating which samples are plotted on the heatmap.
#'          Either a \code{numeric} \code{vector} of indices or a \code{character} \code{vector}
#'          of sample names
#'
#' @param  type \code{character} specifying whether a binary heatmap or a heatmap of continous
#'          deviations is plotted
#'
#' @param threshold \code{numeric} indicates the maximum allowed deviation from means
#'
#' @return
#'
#'
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
#' @import ggplot2 RColorBrewer
#'
#' @export
#'

GC_Heatmap <-function(GcOut,algorithm_step='rt_aligned',substance_subset=NULL,guide='legend',
                      samples_subset=NULL,type="binary",threshold=0.02){

    #########################################
    # A. Select retention times to visualise
    ########################################
    rt_df <- GcOut[[algorithm_step]]
#
#     rt_df <- rt_df[,c(1,3:ncol(rt_df))] # JUst for now, bugfix for rt_extract_heatmap
#     rt_df_org <-rt_df # JUst for now, bugfix for rt_extract_heatmap

    ##########################
    # B Formatting and sorting
    ##########################
    rt_df[,'id'] <- as.character(rt_df[,'id'])
    rt_df <- rt_df[match(as.character(GcOut[["rt_raw"]][,1]),as.character(rt_df[,1])),]

    ##################################################
    # C select a subset of substances, by their position
    ##################################################
    if(!is.null(substance_subset)){
        rt_df <- rt_df[,c(1,substance_subset+1)] # always keep id !, therefore +1
    }
    #############################################
    # D slect a subset of samples, by their names
    #############################################

        if(!is.null(samples_subset)){
            if(is.character(samples_subset)){
        rt_df <- rt_df[rt_df[,1]%in%samples_subset,]
        }else if(is.numeric(samples_subset)){
            rt_df <- rt_df[samples_subset,]
        }
    }
    ############################################
    # E Formatting for ggplot2
    ############################################
    heat_matrix <- reshape2::melt(data = rt_df,id.vars ='id')
    names(heat_matrix) <- c('id','substance','rt')
    heat_matrix[,'substance'] <- as.numeric(as.character(heat_matrix[,'substance']))


    ###########################################
    # F Estimate difference from mean of each row
    ###########################################

    heat_matrix[,'diff'] <- (as.numeric(heat_matrix[,'rt']) - heat_matrix[,'substance'])
    heat_matrix['diff'][heat_matrix['rt']==0] <- 0

    heat_matrix[,'id'] <- ordered( heat_matrix[,'id'], levels = as.factor(rt_df[,'id']))
    heat_matrix[,'substance'] <- ordered( heat_matrix[,'substance'], levels = as.factor(colnames(rt_df)))


    ###################
    # G Binary coding
    ###################

    if(type=="binary"){
        heat_matrix['diff'][abs(heat_matrix['diff'])>threshold] <- 1 # Deviates
        heat_matrix['diff'][abs(heat_matrix['diff'])<threshold] <- 0 # Is okay
        heat_matrix['diff'][heat_matrix['rt']==0] <- NA
        }

    ###################
    # Plot the Heatmap
    ###################

    if(type=="binary"){
        if(max(heat_matrix['diff'],na.rm = T)==0){ # Zeros indicates no deviations of any substance
        gg <- ggplot(heat_matrix, aes(x=substance, y=id,fill=diff),colour="Blue")
        gg <- gg + geom_tile(color="transparent", size=0.001)
        gg <- gg + scale_fill_gradientn(colours = 'blue',na.value = "white")
        gg <- gg + labs(x=NULL, y=NULL, title=paste("No deviations exceeding a threshold of",as.character(threshold)))
        gg <- gg + guides(fill=FALSE)

        }else if(min(heat_matrix['diff'],na.rm = T)==1){ # Really bad alignment
            gg <- ggplot(heat_matrix, aes(x=substance, y=id,fill=diff),colour="red")
            gg <- gg + geom_tile(color="transparent", size=0.001)
            gg <- gg + scale_fill_gradientn(colours = 'red',na.value = "white")
            gg <- gg + labs(x=NULL, y=NULL, title=paste("Alignment failed, all substances show deviation from threshold of",as.character(threshold)))
            gg <- gg + guides(fill=FALSE)

        }else{ # Should be the general outcome, some samples still deviate
            gg <- ggplot(heat_matrix, aes(x=substance, y=id,fill=diff))
            gg <- gg + geom_tile(color="transparent", size=0.001)
            gg <- gg + scale_fill_continuous(low = "yellow",high = "black",breaks=c(0,1),na.value = "white",
                                             guide = 'legend',name=paste('Deviation\n','>',as.character(threshold)),labels=c('NO','YES'))
            gg <- gg + labs(x=NULL, y=NULL, title=paste("Deviations of retention times at a threshold of",as.character(threshold)))

        }
    }else{
        myPalette <-  colorRampPalette(rev(RColorBrewer::brewer.pal(11, "BrBG"))) # take a colour
        gg <- ggplot(heat_matrix, aes(x=substance, y=id, fill=diff))
        gg <- gg + geom_tile(color="white", size=0.01)
        gg <- gg + scale_fill_gradientn(colours = myPalette(10),guide = guide,name='Deviation')
        gg <- gg + labs(x=NULL, y=NULL, title="Deviation of retention times")

        }
    gg <- gg + theme(plot.title=element_text(hjust=0,size = 16,face = 'bold'))
        gg <- gg + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks=element_blank(),
                     axis.text.y=element_text(size = 8))
    gg <- gg + coord_equal(ncol(rt_df)/nrow(rt_df))

    gg
}
