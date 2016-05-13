#' GC_Heatmap visualises the goodness of chromatogram alignment
#'
#' @description Visualises the deviation of retention of sample from the mean of all other,
#'  to indicate how well chromatograms are aligned. Small deviations indicate a good alignment
#'
#' @param GcOut \code{data.frame} representing a matrix of retention times,
#'          where samples are ordered in rows
#'
#' @param all logical, indicating whether all or just the final retention times
#'          are plotted. Default \code{all=FALSE}
#'
#' @return
#'
#'
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import ggplot2
#'
#' @import RColorBrewer
#'
#'
#' @export
#'



GC_Heatmap <-function(GcOut, all=FALSE){

    rt_df <- GcOut[['rt_aligned']]

    ##########################################
    rt_df[,'id'] <- as.character(rt_df[,'id'])

    # Put everything into one 3-column dataframe
    # ##########################################
    heat_matrix <- reshape2::melt(data = rt_df,id.vars ='id')
    names(heat_matrix) <- c('id','substance','rt')
    heat_matrix[,'substance'] <- as.numeric(as.character(heat_matrix[,'substance']))

    ###########################################
    # Estimate difference from mean of each row
    ###########################################

    heat_matrix[,'diff'] <- (as.numeric(heat_matrix[,'rt']) - heat_matrix[,'substance'])
    heat_matrix['diff'][heat_matrix['rt']==0] <- 0


    # Default final alignment
    # optional several plots, which a accessible through click like plot.glm
    # Not included yet

    # format for ggplot, three columns, sample, substance and rt
    # ##########################################################

    steps <- seq(from=1,to=nrow(heat_matrix),by=ncol(rt_df)-1)

    # Orderring of levels
    # ####################
    heat_matrix[,'id'] <-
        ordered( heat_matrix[,'id'], levels = as.factor(rt_df[,'id']))
    suppressWarnings(
    heat_matrix[,'substance'] <-
        ordered( heat_matrix[,'substance'], levels = as.factor(colnames(rt_df)))
    )

    myPalette <- colorRampPalette(rev(brewer.pal(11, "BrBG")))

    # Plot the Heatmap
    # #################
    gg <- ggplot(heat_matrix, aes(x=substance, y=id, fill=diff))
    gg <- gg + geom_tile(color="white", size=0.01)
    gg <- gg + scale_fill_gradientn(colours = myPalette(10),limits=c(-0.4,0.4),
                                    guide = 'colourbar')
    gg <- gg + labs(x=NULL, y=NULL, title="Variation of retention times among samples")
    gg <- gg + theme(plot.title=element_text(hjust=0,size = 22,face = 'bold'))
    # gg <- gg + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    gg <- gg + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks=element_blank(),
                     axis.text.y=element_text(size = 10))
    gg <- gg + coord_equal(10)

    gg




}
# for (i in 1:5){
#     plot(1:i)
#     locator(1)
# }
