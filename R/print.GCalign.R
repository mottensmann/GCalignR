#' print.GCalign: Gives an overview about an executed alignment
#'
#' @details
#' for function arguments that have not been entered the default are printer, but
#' these are update automatically, if these are changed in the main funtion!
#'
#' @param x
#' Object of class \code{GCalign}
#'
#' @param ...
#' actually not useful here, but required to pass the check (S3 method of print)
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#' @export
    print.GCalign <- function(x,...){
# consistent with demands of S3 method, but actually additional arguments
# in addition to x are not executable! Consider revising
cat("Data\n")
    cat(paste("\t",x$call$data,"\n"))
    cat(paste("\t",x$summary$No_of_samples,"samples\n"))
    cat(paste("\t","Reference:",x$call$reference,"\n"))
    if(!is.null(x$call$blanks)){
        temp <- x$call$blanks
        cat("\t","Blanks: ")
        cat(temp,sep = ",")
    }
cat("\nArguments\n")
if(!is.null(x$call$max_linear_shift)){ #check if defined
    cat(paste("\t","Window for linear shifts:",x$call$max_linear_shift,"\n"))
}else{
    cat(paste("\t","Window for linear Shifts:",0.05,"\n"))
}

if(!is.null(x$call$diff_peak2mean)){ #check if defined
    cat(paste("\t","Maximum peak-to-mean distance:",x$call$diff_peak2mean,"\n"))
}else{
    cat(paste("\t","Maximum peak-to-mean distance:",0.02,"\n"))
}

if(!is.null(x$call$diff_peak2peak)){ #check if defined
    cat(paste("\t","Maximum peak-peak distance:",x$summary$diff_peak2peak,"\n"))
}else{
    cat(paste("\t","Maximum peak-peak distance:",0.03,"\n"))
}

cat("\nOutput")
    cat(paste("\t", "Number of distinct peaks:",x$summary$No_Peaks_aligned),"\n")
    if(!is.null(x$call$delete_single_peak)){ #check if defined
        if(x$call$delete_single_peak){
    cat(paste("\t", "Number of singular peaks deleted:",x$summary$singular_peaks,"\n"))
        }
        }

    cat(paste("\t","Range of retention times:",round(min(x$aligned[[x$call$rt_col_name]]$mean_RT),2),
"-",round(max(x$aligned[[x$call$rt_col_name]]$mean_RT),2)),"\n")

    cat(paste("\t","Mean relative Variance in retention time per peak before alignment:",round(x$summary$variance_before$average,4),"\n"))
    cat(paste("\t","Mean relative Variance in retention time per peak after alignment:",round(x$summary$variance_aligned$average,4),"\n"))
}
