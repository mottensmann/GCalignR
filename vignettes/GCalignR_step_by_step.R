## ---- echo = FALSE-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", cache = FALSE,
    fig.width = 6, fig.height = 6) # warning = FALSE

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("mastoffel/GCalignR", build_vignettes = TRUE)

## ------------------------------------------------------------------------
library("GCalignR") 

## ---- eval = FALSE-------------------------------------------------------
#  ?GCalignR

## ---- out.width = 650, fig.retina = NULL, echo = FALSE-------------------
knitr::include_graphics("example.jpg") 

## ---- out.width = 650, fig.retina = NULL, echo = FALSE-------------------
knitr::include_graphics("example2.jpg") 

## ------------------------------------------------------------------------
data("peak_data")
length(peak_data) # number of individuals, i.e. number of list elements
names(peak_data) # names of individuals, i.e. names of list elements 
head(peak_data[[1]]) # column names and data, i.e. one data.frame of list element 

## ------------------------------------------------------------------------
check_input(data = peak_data,list_peaks = F) # If list_peaks = T, a histogram of peaks is plotted 

## ---- eval = FALSE-------------------------------------------------------
#  peak_data <- peak_data[1:4] # subset for speed reasons
#  peak_data_aligned <- align_chromatograms(data = peak_data, # input data
#      conc_col_name = "area", # peak abundance variable
#      rt_col_name = "time", # retention time
#      rt_cutoff_low = 5, # cut peaks with retention times below 5 Minutes
#      rt_cutoff_high = 45, # cut peaks with retention times above 45 Minutes
#      reference = NULL, # Reference will be choosen automatically
#      max_linear_shift = 0.05, # maximum linear shift of chromatograms
#      max_diff_peak2mean = 0.03, # maximum distance of a peak to the mean
#      min_diff_peak2peak = 0.03, # maximum distance between the mean of two peaks
#      blanks = "C2", # no blanks. Specify blanks by names (e.g. c("blank1", "blank2"))
#      delete_single_peak = TRUE, # delete peaks that are present in just one sample
#      write_output = NULL) # add c("time","area") to write data frames to .txt file

## ---- eval = FALSE-------------------------------------------------------
#  peak_data_aligned$aligned$time # to access the aligned retention times
#  peak_data_aligned$aligned$area # to access the aligned area data

## ------------------------------------------------------------------------
data("aligned_peak_data") # the package includes the already aligned data set  

## ----message=FALSE,fig.width=8,fig.height=8------------------------------
gc_heatmap(aligned_peak_data,threshold = 0.05) # By default a threshold of 0.05 is used to mark deviations

