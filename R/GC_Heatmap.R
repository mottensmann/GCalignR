#' GC_Heatmap visualises the goodness of chromatogram alignment
#'
#' @description Visualises the deviation of retention of sample from the mean of all other,
#'  to indicate how well chromatograms are aligned. Small deviations indicate a good alignment
#'
#' @param GcOut \code{data.frame} representing a matrix of retention times,
#'          where samples are ordered in rows
#'
#' @param step character, indicating which step of the algorithm to plot. Either
#'          {rt_raw}, {rt_linear} or {rt_aligned}
#'
#' @param substance_subset vector containing indices of substances (i.e. rows) to plot
#'          By default {NULL} indicating all substances are plotted
#'
#' @param guide character, indicating type of colourbar as discrete (i.e 'legend')
#'          or gradient (i.e 'colourbar)
#'
#' @param  limits vector, allows to finetune the scaling of the colourbar
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



GC_Heatmap <-function(GcOut,step='rt_aligned',substance_subset=NULL,guide='colourbar',
                      limits=c(-0.05,0.05)){
    rt_df <- GcOut[[step]]

    if(!is.null(substance_subset)){
        substance_subset <- substance_subset+1 # because first holds ids

        rt_df <- rt_df[,c(1,substance_subset)] # always keep id
    }

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
    gg <- gg + scale_fill_gradientn(colours = myPalette(10),limits=limits,
                                    guide = guide)
    gg <- gg + labs(x=NULL, y=NULL, title="Variation of retention times among samples")
    gg <- gg + theme(plot.title=element_text(hjust=0,size = 22,face = 'bold'))
    # gg <- gg + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    gg <- gg + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks=element_blank(),
                     axis.text.y=element_text(size = 6))
    gg <- gg + coord_equal(ratio = ncol(rt_df)/nrow(rt_df))

    gg




}
# for (i in 1:5){
#     plot(1:i)
#     locator(1)
# }
