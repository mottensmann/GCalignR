#' Visualise deviation of retention times from the mean per row
#' gc_heatmap visualises the current variation of retention times among all samples
#' within the same row. Therefore the absolute deviation from the mean of the row is calculated
#' for each sample
#'
#'
#' @param rts_matrix Matrix of Retention times, where samples are sorted in rows
#'          and retention times in columns. Absence of substances in a sample are
#'          coded by NA. These values are ignored in calculating the row means
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


GC_Heatmap <-function(rts_matrix){

    # Default final alignment
    # optional several plots, which a accessible through click like plot.glm
    # Not included yet

    # format for ggplot, three columns, sample, substance and rt
    # ##########################################################

    heat_matrix <- matrix(data = NA,nrow = nrow(rts_matrix)*ncol(rts_matrix),ncol = 3)
    steps <- seq(from=1,to=nrow(rts_matrix)*ncol(rts_matrix),by=ncol(rts_matrix))

    # Obtain average retention times of substances
    # ############################################
    comp_means <- colMeans(rts_matrix,na.rm=T)

    # Put everything into one 3-column dataframe
    # ##########################################

    for (i in 1:nrow(rts_matrix)){
heat_matrix[steps[i]:(steps[i]+(ncol(rts_matrix)-1)),1] <- rownames(rts_matrix)[i]
heat_matrix[steps[i]:(steps[i]+(ncol(rts_matrix)-1)),2] <- (rts_matrix[i,]-comp_means) # before abs()
heat_matrix[steps[i]:(steps[i]+(ncol(rts_matrix)-1)),3] <- colnames(rts_matrix)
    }
    heat_matrix <- data.frame(ID = heat_matrix[,1],
                              val = as.numeric(as.character(heat_matrix[,2])),
                              RT = heat_matrix[,3])
    # Zero deviation, if sampe is empty at the specific substance
    #############################################################
    heat_matrix[is.na(heat_matrix)] <- 0

    # Orderring of levels
    # ####################
    heat_matrix[,3] <-
    ordered( heat_matrix[,3], levels = as.factor(colnames(rts_matrix)))
    heat_matrix[,1] <-
        ordered( heat_matrix[,1], levels = as.factor(rownames(rts_matrix)))

    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

    # Plot the Heatmap
    # #################
    gg <- ggplot(heat_matrix, aes(x=RT, y=ID, fill=val))
    gg <- gg + geom_tile(color="white", size=0.01)
    gg <- gg + scale_fill_gradientn(colours = myPalette(3),limits=c(-0.4,0.4))
    gg <- gg + labs(x=NULL, y=NULL, title="Variation of Retention times per Substance among Chromatograms")
    gg <- gg + theme(plot.title=element_text(hjust=0))
    # gg <- gg + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    gg <- gg + theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks=element_blank())
    gg <- gg + coord_equal(10)
    gg

}



