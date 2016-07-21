#' Plot Diagonstics for an Gcalign Object
#'
#' @description
#' Two plots are currently available: One plot visualises the distribution of linear shifts
#' that were applied to align chromatograms to a reference before aligning individual peaks.
#' A second plot illustrates the remaining variation of retention times on the level of individual
#' peaks by plotting the distribution of retention time ranges.
#'
#' @usage
#' ## S3 method for class "GCalign"
#' plot(object)
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
    lin_shift_table <- function(x){
        # 1.) Get the search window used in the function call and setup a table with all steps
        # 2.) Count the applied steps in the Alignment process
        # 3.) calculate a frequency table
        x <- as.data.frame(object[["Logfile"]][["Call"]]) # function call
        x <- as.numeric(as.character(x[["max_linear_shift"]])) # get the xdow for linear shifts
        x <- seq(-x,x,0.01) # all linear steps investigated

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

    object_to_matrix <- function(object,step="aligned",rt_col="RT"){
        L <- length(object[[step]][[rt_col]])-1
        rt_mat <- matrix(data = NA,nrow = L,ncol = length(object[["aligned"]][[rt_col]][[1]]))
        for(i in 1:L){
            rt_mat[i,] <- object[[step]][[rt_col]][[i+1]]
        }
        return(rt_mat)
    }

    bw <- function(b,x){b/bw.nrd0(x)} # used as a helper to smooth the gaussian fit

     multiplot <- function(..., plotlist=NULL, cols) { # cookbook.R
        require(grid)
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        numPlots = length(plots)
        # Make the panel
        plotCols = cols # Number of columns of plots
        plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
        vplayout <- function(x, y)
            viewport(layout.pos.row = x, layout.pos.col = y)
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            curRow = ceiling(i/plotCols)
            curCol = (i-1) %% plotCols + 1
            print(plots[[i]], vp = vplayout(curRow, curCol ))
        }
    }
    # -----------------------------------------------------------------------------------------
df <- lin_shift_table(object) # steps of linear shifts and their frequency
long <- nrow(df) # x pos. for annotation in ggplot
lat <- max(df["count"]) # y pos.

LinShift <- ggplot2::ggplot(data = df,aes(shift,count)) +
    geom_bar(stat = "identity",fill="blue") +
    labs(title ="Linear Adjustments",
                x = "Shift Size",
                y = "Frequency") +
    theme_bw()+theme(
    plot.title=element_text(face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
    panel.grid.major = element_blank()) +
    geom_segment(aes(x=1,y=lat+0.01,xend=long,yend=lat+0.01),arrow = arrow(length = unit(0.02,"npc")),color="red",size=0.8)+
    geom_segment(aes(x=long,y=lat+0.01,xend=1,yend=lat+0.01),arrow = arrow(length = unit(0.02,"npc")),color="red",size=0.8)+
    annotate("text",y=lat+0.02,x=round(nrow(df)/2),label=paste0("Window"),size=5)
#**********************************************************************************
aligned <- MinMax(object[["heatmap_input"]][["aligned_rt"]][,-1]) # Range of RTs aligned
aligned <- as.data.frame(aligned["range"]) # Formatting

RT_Range <-    ggplot() +
    geom_histogram(aes(x=range,y=..density..),data = aligned,binwidth = 0.01,colour="black",fill="white")+
    geom_density(aes(x=range),data = aligned,adjust=bw(0.006,aligned[["range"]]),alpha=0.6,fill="Blue")+
    labs(title ="Variation of Retention Times per Peak",
         x = "Range",
         y = "Frequency")+
    theme_bw()+theme(
        plot.title=element_text(face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        axis.ticks.y = element_line(size = 0.5, colour = "grey40"),
        axis.ticks.x = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank())


multiplot(LinShift,RT_Range,cols = 2)
}


