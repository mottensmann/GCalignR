# RangePlot


range_hist <- function(GCalign,Aligner=NULL){
    if(!class(GCalign)=="GCalign")stop("Data is not of class 'GCalign'")
    aligned <- MinMax(GCalign$heatmap_input$aligned_rt[,-1])
    lintrans <- MinMax(GCalign$heatmap_input$linear_shifted_rt[,-1])
    int <- MinMax(GCalign$heatmap_input$initial_rt[,-1])
    out <- matrix(NA,nrow = nrow(aligned),ncol=2) #3
    colnames(out) <- c("Linear Adjustment","Input") # c("aligned","linshift","input")
    out <- as.data.frame(out)
    #out[1:nrow(aligned),1] <- aligned$range
    out[1:nrow(lintrans),1] <- lintrans$range #2
    out[1:nrow(int),2] <- int$range#3
    df.m <- reshape2::melt(out)

    if(!is.null(Aligner)){
    p <- ggplot(aligned)
    p <- p + stat_density(aes(x=range,y=..scaled..),fill="Blue")
    } else{
    p <- ggplot(df.m) + stat_density(aes(x=value, y=..scaled..,col=variable),alpha=0.4, geom="line",position = "dodge",size=2)
    p + stat_density(aes(fill=variable))
    }
    p <- p + labs(title ="Variation of Retention Times",
                    x = "Range of Retention Times",
                    y = "Frequency")
    p <- p + theme_bw()+theme(
        plot.title=element_text(face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank())
    p
    return(p)
}

MinMax <- function(rt_mat = aligned){ # Estimate the range of retention times per substance, they should be no overlapp
    temp <- matrix(NA,1,2)
    colnames(temp) <- c("range","row")
    data <- temp[0,]
    for(i in 1:ncol(rt_mat)){
        data<-rbind(data,cbind(abs(diff(range(rt_mat[,i][rt_mat[,i]>0],na.rm = TRUE))),i))
    }
    return(as.data.frame(data))
}

PeakOverlap <- function(rt_min_max,threshold=0.02){
    rt_min_max["width"] <- rt_min_max["max"] - rt_min_max["min"]
    rt_min_max["distance"] <- 0
    for(n in 2:nrow(rt_min_max)){
        rt_min_max[n,"distance"] <- rt_min_max[n,"min"] - rt_min_max[n-1,"max"]
    }
    rt_min_max["distance"][ rt_min_max["distance"]>threshold] <- threshold
    rt_min_max["offset"] <- 0
    for(n in 2:nrow(rt_min_max) ){
        rt_min_max[n,"offset"] <- rt_min_max[n,"distance"] + rt_min_max[n-1,"width"]
    }
    rt_min_max["cumul"] <- cumsum(rt_min_max["distance"]+rt_min_max["offset"])
    return(rt_min_max)

}

GCalign_to_matrix <- function(GCalign,step="aligned",rt_col="RT"){
    L <- length(GCalign[[step]][[rt_col]])-1
    rt_mat <- matrix(data = NA,nrow = L,ncol = length(GCalign[["aligned"]][[rt_col]][[1]]))
    for(i in 1:L){
        rt_mat[i,] <- GCalign[[step]][[rt_col]][[i+1]]
    }
    return(rt_mat)
}
