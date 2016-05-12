#' Visualise deviation of retention times from the mean per row
#' gc_heatmap visualises the current variation of retention times among all samples
#' within the same row. Therefore the absolute deviation from the mean of the row is calculated
#' for each sample
#'
#'
#' @param rt_df Dataframe
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



GC_Heatmap <-function(rt_df){

    rt_df <- x
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
    heat_matrix[,'substance'] <-
        ordered( heat_matrix[,'substance'], levels = as.factor(colnames(rt_df)))

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





#
# GC_Heatmap <-function(rts_df){
#
#
#     # 1. NAs if 0
#     # 2. reshape
#     # 3. differences, 3 decimals
#     # 4. plot
#
#     rts_matrix <- test[['RT']]
#     #library(reshape2)
#
#     rts_matrix <- as.data.frame(t(rts_matrix))
#     rts_matrix <- cbind(rownames(rts_matrix),rts_matrix)
#     colnames(rts_matrix)[2:ncol(rts_matrix)] <- as.character(round(as.numeric(rts_matrix[1,2:ncol(rts_matrix)]),2))
#     colnames(rts_matrix)[1] <- 'id'
#     rts_matrix[,1] <- as.character(rts_matrix[,1])
#
#     heat_matrix <- reshape::melt(data = rts_matrix[2:nrow(rts_matrix),],id.vars ='id')
#
#
#     # Default final alignment
#     # optional several plots, which a accessible through click like plot.glm
#     # Not included yet
#
#     # format for ggplot, three columns, sample, substance and rt
#     # ##########################################################
#
#     steps <- seq(from=1,to=nrow(rts_matrix)*ncol(rts_matrix),by=nrow(rts_matrix))
#
#     # Obtain average retention times of substances
#     # ############################################
#     rts_matrix[is.na(rts_matrix)] <- NA
#     comp_means <- rowMeans(rts_matrix,na.rm=T)
#     rts_matrix[is.na(rts_matrix)] <- NA
#
#     # Put everything into one 3-column dataframe
#     # ##########################################
#
#     for (i in 1:ncol(rts_matrix)){
# heat_matrix[steps[i]:(steps[i]+(nrow(rts_matrix)-1)),1] <- colnames(rts_matrix)[i]
#
# heat_matrix[steps[i]:(steps[i]+(nrow(rts_matrix)-1)),2] <- (rts_matrix[,i]-comp_means) # before abs()
# heat_matrix[steps[i]:(steps[i]+(nrow(rts_matrix)-1)),3] <- as.character(comp_means)
#     }
#     heat_matrix <- data.frame(ID = heat_matrix[,1],
#                               val = as.numeric(as.character(heat_matrix[,2])),
#                               RT = heat_matrix[,3])
#     # Zero deviation, if sampe is empty at the specific substance
#     #############################################################
#     heat_matrix[is.na(heat_matrix)] <- 0
#
#     # Orderring of levels
#     # ####################
#     heat_matrix[,3] <-
#     ordered( heat_matrix[,3], levels = as.factor(colnames(rts_matrix)))
#     heat_matrix[,1] <-
#         ordered( heat_matrix[,1], levels = as.factor(rownames(rts_matrix)))
#
#     myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#
#     # Plot the Heatmap
#     # #################
#     gg <- ggplot(heat_matrix, aes(x=RT, y=ID, fill=val))
#     gg <- gg + geom_tile(color="white", size=0.01)
#     gg <- gg + scale_fill_gradientn(colours = myPalette(3),limits=c(-0.4,0.4))
#     gg <- gg + labs(x=NULL, y=NULL, title="Variation of Retention times per Substance among Chromatograms")
#     gg <- gg + theme(plot.title=element_text(hjust=0))
#     # gg <- gg + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#     gg <- gg + theme(axis.title.x=element_blank(),
#                      axis.text.x=element_blank(),
#                      axis.ticks=element_blank())
#     gg <- gg + coord_equal(10)
#     gg
#
# }



