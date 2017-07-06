## ---- echo = FALSE-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, comment = ">", cache = FALSE,
    fig.width = 8, fig.height = 6, fig.align = "center") 

## ---- results='hide', echo=FALSE-----------------------------------------
library(GCalignR)
library(ggplot2)

## ---- fig.cap="Figure 1. A Chromatogram plots an intensity signal over the course of a separation run.", echo=F----
set.seed(123)
# create one chromatogram with specified peaks
df1 <- GCalignR:::simple_chroma(peaks = c(5.01,10.02,13.10,20.22,24.57), N = 1)
df1 <- subset(df1, x > 3.5 & x < 27)
# plot the chromatogram
chroma <- ggplot(data = df1, aes(x,y, fill = sample)) + geom_line(size = 1) + theme_classic() + xlab("Retention time ") + ylab("Intensity") + scale_x_continuous(breaks = seq(4,26,1),expand = c(0,0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
chroma

## ---- echo=F, fig.cap="Figure 2. Chromatogram with integrated peaks"-----
# Using an internal function, peaks are detected by searching for global maxima
peaks <- find_peaks(df1) 
# create the plot
chroma + geom_linerange(data = peaks, aes(x = x, ymin = 0, ymax = y), linetype = "dashed", col = "black") + annotate("text", x = peaks[["x"]], y = peaks[["y"]] + 0.1, label = as.character(round(peaks[["y"]],2)), angle = 0) + geom_area(fill = "blue", alpha = 0.4) + theme(legend.position = "none")

## ---- echo=FALSE---------------------------------------------------------
# create a data frame that depicts a peak list
df <- data.frame(row.names = c("Peak 1", "Peak 2", "Peak 3", "Peak 4", "Peak 5"), time = peaks[["x"]], height = peaks[["y"]])
# print the table
knitr::kable(df,digits = 2)

## ---- results='hide', echo = FALSE, fig.cap="Figure 3. Overlay of Chromatograms from four samples"----
set.seed(123)
peak_list <- sample(x = seq(from = 1, to = 26, by = 4), size = 6, replace = F)
df <- GCalignR:::simple_chroma(peaks = peak_list, N = 4)
# draw chromatograms and display peaks
chroma <- ggplot(data = df, aes(x,y, col = sample)) + geom_line(size = 1) + theme_classic() + xlab("Retention time ") + ylab("Intensity") + scale_x_continuous(breaks = seq(0,30,5),expand = c(0,0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom") + scale_color_brewer(palette = "Dark2") + guides(col = guide_legend(ncol = 4, title = NULL))

# peaks are in this case simply the local maxima for each sample
peaks <- find_peaks(df)
chroma <- chroma + geom_linerange(data = peaks, aes(x = x, ymin = y, ymax = y + 0.1), linetype = "solid", col = "black") + annotate("text", x = peaks[["x"]], y = peaks[["y"]] + 0.2, label = as.character(round(peaks[["x"]],2)), angle = 90)
print(chroma)

## ---- echo = FALSE, results='hide', eval=FALSE---------------------------
#  ## these lines create a input file which is distributed with the package
#  # sink("ChromSimul.txt",append = FALSE)
#  ## write sample identifier
#  #cat(levels(peaks[["sample"]]),sep = "\t")
#  ## write variables
#  #cat(c("\nrt","height\n"),sep = "\t")
#  ## merge data horizontally
#  
#  #dat_mat <- numeric()
#  
#  #for (i in levels(peaks[["sample"]])) {
#   #   temp <- as.matrix(peaks[,c("x","y")][peaks[["sample"]] == i,])
#    #  add <-  max(summary(peaks[["sample"]])) - nrow(temp)
#     # temp <- rbind(temp, matrix(data = 0,nrow = add, ncol = 2))
#      #  dat_mat <- cbind(dat_mat, temp)
#  #}
#  #write.table(dat_mat, row.names = F, col.names = F, sep = "\t")
#  #sink()

## ---- echo = FALSE-------------------------------------------------------
## sample identifiers
cat(levels(peaks[["sample"]]),sep = "\t") 
## variable names
cat(c("\nrt","height\n"),sep = "\t") 

## empty matrix to fill with the data
dat_mat <- numeric()
for (i in levels(peaks[["sample"]])) {
    temp <- as.matrix(peaks[,c("x","y")][peaks[["sample"]] == i,])
    add <-  max(summary(peaks[["sample"]])) - nrow(temp)
    temp <- rbind(temp, matrix(data = 0,nrow = add, ncol = 2))
    dat_mat <- cbind(dat_mat, temp)
} 
## output the formatted matrix
write.table(round(dat_mat,2), row.names = F, col.names = F, sep = "\t")

## ---- echo = FALSE, fig.cap="Figure 4. Corrected linear drift between Reference and Focal sample at a shift of -1"----
# create two chromatograms
df1 <- data.frame(simple_chroma(peaks = c(5,10,16,20,24),
                                N = 1,
                                Names = "Reference",
                                sd = c(0.30, 0.35, 0.24, 0.25, 0.23)),
                                ym1 = NA, ym2 = NA, yp1 = NA, yp2 = NA)

x <- c(6,11.2,20.95)
df2 <- data.frame(simple_chroma(peaks = x,
                                N = 1,
                                Names = "Focal sample",
                                sd = c(0.22, 0.28, 0.27)),
                  ym1 = data.frame(simple_chroma(peaks = x - 1,
                                                 N = 1,
                                                 Names = "Focal sample",
                                                 sd = c(0.22, 0.28, 0.27)),
                                   type = "Focal sample")[["y"]],
                  ym2 = data.frame(simple_chroma(peaks = x - 2,
                                                 N = 1,
                                                 Names = "Focal sample",
                                                 sd = c(0.22, 0.28, 0.27)),
                                   type = "Focal sample")[["y"]],
                  yp1 = data.frame(simple_chroma(peaks = x + 1,
                                                 N = 1,
                                                 Names = "Focal sample",
                                                 sd = c(0.22, 0.28, 0.27)),
                                   type = "Focal sample")[["y"]],
                  yp2 = data.frame(simple_chroma(peaks = x + 2,
                                                 N = 1,
                                                 Names = "Focal sample",
                                                 sd = c(0.22, 0.28, 0.27)),
                                   type = "Focal sample")[["y"]])
df <- rbind(df1,df2)

# reduce range of x
df <- subset(df, x > 3.5 & x < 27)
# find peaks
peaks <- find_peaks(df[,1:3])
peaks2 <- peaks
peaks2[["x"]][6:8] <- peaks2[["x"]][c(1:2,4)]

tx <- peaks[["x"]][peaks[["sample"]] == "Focal sample"]
tex <- data.frame(x = c(tx - 2, tx - 1, tx, tx + 1, tx + 2), y = rep(peaks[["y"]][peaks[["sample"]] == "Focal sample"],5), z = rep(-2:2, each = 3), sample = "Focal sample")


chroma <- ggplot(data = df, aes(x,y, fill = sample)) + geom_line(size = 1.2, colour = "black") +
    theme_classic(base_size = 12) + xlab("Retention time ") + ylab("") +
    scale_x_continuous(breaks = seq(4,26,1),expand = c(0,0)) +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
          strip.background = element_rect(colour = "black", fill = "#CCCCFF")) +
    facet_wrap(~sample, ncol = 1) +
    scale_fill_manual(values = c("#1B9E77","#F16913")) +
    theme(legend.position = "none", axis.text.x = element_blank()) +
    geom_area(alpha = 0.5) +
    geom_line(aes(x = x, y = ym1), colour = "#FD8D3C", size = 1.2, linetype = "solid") +
    geom_line(aes(x = x, y = ym2), colour = "#FDAE6B", size = 0.8, linetype = "solid") +
    geom_line(aes(x = x, y = yp1), colour = "#D94801", size = 0.8, linetype = "solid") +
    geom_line(aes(x = x, y = yp2), colour = "#A63603", size = 0.8, linetype = "solid") +
    geom_text(data = tex, aes(x = x, y = y + 0.1), label = as.character(tex[["z"]])) +
    geom_rect(data = peaks2, aes(xmin = x - 0.1, xmax = x + 0.1, ymin = 0, ymax = y), linetype = "solid", colour = "black", fill = "grey50", size = 0.8, alpha = 0.2)
print(chroma)

## ---- echo = FALSE, fig.cap="Figure 4. Corrected linear drift between reference and sample", eval=FALSE----
#  # REPLACED BY EXAMPLE ABOVE
#  # graphical representation of the procedure that is applied during the correction of linear drift.
#  set.seed(1533)
#  # create two samples
#  df <- rbind(simple_chroma(peaks = c(10,15,22), N = 1, Names = "Reference"), simple_chroma(peaks = c(12,16.89,24), N = 1, Names = "Sample"))
#  peaks <- find_peaks(df)
#  peaks2 <- peaks
#  # Adjust retention times to shift the sample chromatogram
#  peaks2[["x"]] <- peaks2[["x"]] - 2
#  peaks2[["x"]][peaks2[["x"]] < 0] <- 0
#  
#  # create a data frame for plotting
#  df2 <- data.frame(x = rep(seq(0,30,length = 10000),2), y = c(GCalignR:::peaks2chroma(data = peaks, sample = "Reference"), GCalignR:::peaks2chroma(data = peaks, sample = "Sample")), sample = rep(c("Reference", "Sample A"), each = 10000), y2 = c(rep(0, 10000),GCalignR:::peaks2chroma(data = peaks2, sample = "Sample")))
#  
#  #subset on x
#  df2 <- subset(df2, df2$x > (min(peaks[["x"]]) - 1) & df$x < (max(peaks[["x"]]) + 1))
#  
#  # plot
#  plot <- ggplot(df2, aes(x = x, y = y, col = sample, fill = sample)) +
#      geom_line(size = 1.2) +
#      theme_classic() + xlab("Retention time ") + ylab("Intensity") +
#      scale_x_continuous(expand = c(0,0)) +
#      theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") +
#      geom_line(aes(x = x,y = y2), linetype = "dotted", size = 1) +
#      scale_color_manual(values = c("black","black"), guide = FALSE) +
#      geom_area() +
#      scale_fill_manual(values = c("darkorange","blue"), name = "") +
#      geom_area(aes(x = x,y = y2), alpha = 0.6) +
#      geom_segment(aes(x = peaks[peaks[["sample"]] == "Sample",][2,1], xend = peaks2[peaks2[["sample"]] == "Sample",][2,1], y = find_peaks(df2)[5,2] + 0.05, yend = find_peaks(df2)[5,2] + 0.05),size = 1.2, colour = "black", arrow = arrow(length = unit(x = 0.14, units = "cm"))) +
#      annotate("text", x = find_peaks(df2)[5,1] - 1, y = find_peaks(df2)[5,2] + 0.25, label = "Linear shift", angle = 0) + xlab(label = "") + ylab(label = "")
#  print(plot)

## ---- echo=FALSE, fig.cap="Figure 5. Alignment of individual peaks based on retention time matrices. Colours represent substances, black rectangles highlight causes of manipulations.",out.width = "750px"----
knitr::include_graphics("align_peaks.png",dpi = 300)

## ---- fig.cap="Figure 5. Chromatographic representation of the dataset prior to alignment"----
## path to the data
path <- system.file("extdata", "simulated_peak_data.txt", package = "GCalignR")
## draw chromatograms
x <- draw_chromatogram(data = path, rt_col_name = "rt", show_rt = T, show_num = F, plot = F)
x[["ggplot"]] + geom_line(size = 1.2) + theme(axis.ticks.x = element_blank()) + ggplot2::scale_color_brewer(palette = "Dark2")

## ---- eval=T, results="hide"---------------------------------------------
aligned <- align_chromatograms(data = path,
                               rt_col_name = "rt",
                               max_linear_shift = 2,
                               max_diff_peak2mean = 0.02,
                               min_diff_peak2peak = 1,
                               reference = "A2")

## ------------------------------------------------------------------------
print(aligned[["Logfile"]][["LinearShift"]])

## ---- results="hide"-----------------------------------------------------
x <- draw_chromatogram(data = aligned, rt_col_name = "rt", step = "lin_aligned", show_rt = F, show_num = F, plot = F)
x[["ggplot"]] + ggplot2::scale_color_brewer(palette = "Dark2")

## ---- results="hide"-----------------------------------------------------
x <- draw_chromatogram(data = aligned, rt_col_name = "rt", step = "fully_aligned", show_num = T, plot = F)
x[["ggplot"]] + ggplot2::scale_color_brewer(palette = "Dark2")

## ------------------------------------------------------------------------
## for using ggplot2::facet_wrap we need to get rid of the annotations
x <- draw_chromatogram(data = aligned, rt_col_name = "rt", step = "fully_aligned", show_num = F, plot = F)
x[["ggplot"]] + ggplot2::facet_wrap(~sample, ncol = 1) + ggplot2::scale_color_brewer(palette = "Dark2")

