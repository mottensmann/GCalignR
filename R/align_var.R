#' Quantify goodness of alignment by relative variances of peak retention times
#'
#' @description
#' \code{align_var()} calculates range and mean of the coefficient of variation for retention
#' times of peaks that are currently within the same row in data.frames of \code{gc_peak_list}
#' (i.e. classified as same substance).
#'
#' @inheritParams matrix_append
#'
#' @inheritParams align_chromatograms
#'
#' @return
#' a list of two numeric vectors
#' \item{range}{the minimum and maximum relative variation of peak retention times}
#' \item{average}{the mean relative variation of peak retention times}
#'
#' @keywords internal
#' @export
#'

# align_var <- function(gc_peak_list,rt_col_name){
#
# variance_in_peaks <- function(gc_peak_list, sample_indices, peak, rt_col_name){
#     rt <- unlist(lapply(gc_peak_list[sample_indices], function(x) x[peak, rt_col_name]))
#     var_rt <- stats::sd(rt[!(rt == 0)], na.rm = TRUE)
#     rel_var_rt <- var_rt/mean(rt,na.rm = TRUE)
#     rel_var_rt
# }
#     peak <- nrow(gc_peak_list[[1]]) # number of peaks
#     out <- unlist(lapply(1:peak,
#     function(x) variance_in_peaks(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
#
#     output <- list(range=range(out,na.rm = T),average=mean(out,na.rm=T))
#
#     return(output)
# }

align_var <- function(gc_peak_list,rt_col_name){
# Calculates the range of retention times for each peak, estimate is the width
# computated as the distance between min and max values
width_of_peaks <- function(gc_peak_list, sample_indices, peak, rt_col_name){
    rt <- unlist(lapply(gc_peak_list[sample_indices], function(x) x[peak, rt_col_name])) # all rts per peak
    if(any(!is.na(rt))){ # Check if all are empty
    min2max <- range(rt[!(rt == 0)], na.rm = TRUE) # range of rts, width per peak
    width <- abs(diff(min2max)) # If zero, no deviation, or just one sample!
    } else{
    width <- NA
    }
    return(width)
    }
peak <- nrow(gc_peak_list[[1]]) # number of peaks
out <- unlist(lapply(1:peak,
 function(x) width_of_peaks(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
output <- list(range=round(range(out,na.rm = T),2),average=round(mean(out,na.rm=T),2),std=round(sd(out,na.rm = T),2))
return(output)
}


