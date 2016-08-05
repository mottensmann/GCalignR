#' Plot Diagonstics for an Gcalign Object
#'
#' @description
#' Three plots are currently available: One plot visualises the distribution of linear shifts
#' that were applied to align chromatograms to a reference before aligning individual peaks.
#' A second plot illustrates the remaining variation of retention times on the level of individual
#' peaks by plotting the distribution of retention time ranges. The third plots shows a distribution
#' of peak numbers after aligning the chromatograms.
#'
#' @usage
#' ## S3 method for class "GCalign"
#' plot(x,...)
#'
#' @return
#' a ggplot2 figure including three subplots
#'
#' @param x
#' \code{GCalign} object, result of \code{\link{align_chromatograms}}
#'
#' @param ...
#' optional arguments
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @export
#'
plot.GCalign <- function(x,...){

    # Define  internal function
    # -----------------------------------------------------------------------------------------
    lin_shift_table <- function(object){
        # 1.) Get the search window used in the function call and setup a table with all steps
        # 2.) Count the applied steps in the Alignment process
        # 3.) calculate a frequency table
        x <- as.data.frame(object[["Logfile"]][["Call"]]) # function call
        x <- as.numeric(as.character(x[["max_linear_shift"]])) # get the xdow for linear shifts
        x <- seq(-x[[1]],x[[1]],0.01) # all linear steps investigated

        df <- matrix(0,nrow = length(x),ncol = 2) # matrix of steps and counts
        df <- as.data.frame(df,row.names = F)
        names(df) <- c("shift","count")
        df["shift"] <- as.factor(x) # template to be filled up with applied shifts

        x <- as.data.frame(table(object[["Logfile"]][["LinearShift"]]["shift"])) # x Linshift data.frame
        names(x) <- c("shift","count")
        x["count"] <- x["count"] / sum(x["count"])
        index <- which((df[,1] %in% x[,1])) # which steps in df to overwrite
        df[index,"count"] <- x["count"] # combined tables
        df["shift"] <- as.factor(round(as.numeric(as.character(df[["shift"]])),2))
        return(df)
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

    object_to_matrix <- function(x,step="aligned",rt_col="RT"){
        L <- length(x[[step]][[rt_col]])-1
        rt_mat <- matrix(data = NA,nrow = L,ncol = length(x[["aligned"]][[rt_col]][[1]]))
        for(i in 1:L){
            rt_mat[i,] <- x[[step]][[rt_col]][[i+1]]
        }
        return(rt_mat)
    }

    bw <- function(b,x){b/bw.nrd0(x)} # used as a helper to smooth the gaussian fit

         # -----------------------------------------------------------------------------------------
df <- lin_shift_table(x) # steps of linear shifts and their frequency
long <- nrow(df) # x pos. for annotation in ggplot
lat <- max(df["count"]) # y pos.

LinShift <- ggplot2::ggplot(data = df,aes(shift,count)) +
    geom_bar(stat = "identity",fill="navyblue") +
    labs(title ="Linear Transformation of Retention Times",
                x = "Shift",
                y = "Frequency") +
    theme_bw()+theme(
    plot.title=element_text(face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 16,angle=45,hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
    panel.grid.major = element_blank())
    # geom_segment(aes(x=1,y=lat+0.01,xend=long,yend=lat+0.01),arrow = arrow(length = unit(0.02,"npc")),color="red",size=0.8)+
    # geom_segment(aes(x=long,y=lat+0.01,xend=1,yend=lat+0.01),arrow = arrow(length = unit(0.02,"npc")),color="red",size=0.8)+
    # annotate("text",y=lat+0.02,x=round(nrow(df)/2),label=paste0("Window"),size=5)
#**********************************************************************************
aligned <- MinMax(x[["heatmap_input"]][["aligned_rt"]][,-1]) # Range of RTs aligned
aligned <- as.data.frame(aligned["range"]) # Formatting

RT_Range <-    ggplot() +
    #geom_histogram(aes(x=range,y=..density..),data = aligned,binwidth = 0.01,colour="black",fill="gray85")+
    geom_density(aes(x=range),data = aligned,adjust=bw(0.006,aligned[["range"]]),alpha=1,fill="limegreen")+
    labs(title ="Variation of Retention Times per Peak",
         x = "Range",
         y = "Frequency")+
    theme_bw()+theme(
        plot.title=element_text(face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 16,angle = 45,hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank())
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

peaks_final <- ggplot2::ggplot(peak_df,aes(x=ID,y=Peaks)) +
    ggplot2::geom_bar(stat = "identity",fill="navyblue") +
    theme_bw()+theme(
        plot.title=element_text(face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x=element_text(angle=90,vjust = +0.5),
        axis.ticks.x = element_line(size = 1,colour = "black"),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank())+
labs(title ="Posterior Distribution of Peaks",
         x = "",
         y = "Number of Peaks")+
    geom_text(aes(label=Peaks),size=3,colour="gray5", position=position_dodge(width=0.5), vjust=-0.25)


return(gridExtra::grid.arrange(gridExtra::arrangeGrob(LinShift, RT_Range,nrow = 1),peaks_final, nrow = 2))

}


