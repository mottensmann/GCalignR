#' print.GCalign: Gives an overview about an executed alignment
#'
#' @details
#' for function arguments that have not  been entered the default are printer, but
#' these are update automatically, if these are changed in the main funtion !
#'
#' @param object
#' Object of class \code{GCalign}
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#' @export
    print.GCalign <- function(object){

cat("Data\n")
    cat(paste("\t",object$call$data,"\n"))
    cat(paste("\t",object$summary$No_of_samples,"samples\n"))
    cat(paste("\t","Reference:",object$call$reference,"\n"))
    if(!is.null(object$call$blanks)){
        temp <- object$call$blanks
        cat("\t","Blanks: ")
        cat(temp,sep = ",")
    }
cat("\nArguments\n")
if(!is.null(object$call$max_linear_shift)){ #check if defined
    cat(paste("\t","Window for linear shifts:",object$summary$maobject_diff_peak2mean,"\n"))
}else{
    cat(paste("\t","Window for linear Shifts:",0.05,"\n"))
}

if(!is.null(object$call$diff_peak2mean)){ #check if defined
    cat(paste("\t","Favoured peak-to-mean distance:",object$summary$maobject_diff_peak2mean,"\n"))
}else{
    cat(paste("\t","Favoured peak-to-mean distance:",0.02,"\n"))
}

if(!is.null(object$call$diff_peak2peak)){ #check if defined
    cat(paste("\t","Favoured peak-peak distance:",object$summary$maobject_diff_peak2mean,"\n"))
}else{
    cat(paste("\t","Favoured peak-peak distance:",0.03,"\n"))
}



cat("\nOutput")
    cat(paste("\t", "Number of distinct peaks:",object$summary$No_Peaks_aligned),"\n")
    if(!is.null(object$call$delete_single_peak)){ #check if defined
        if(object$call$delete_single_peak){
    cat(paste("\t", "Number of singular peaks deleted:",object$summary$singular_peaks,"\n"))
        }
        }

    cat(paste("\t","Range of retention times:",round(min(object$aligned[[object$call$rt_col_name]]$mean_RT),2),
"-",round(max(object$aligned[[object$call$rt_col_name]]$mean_RT),2)),"\n")

    cat(paste("\t","Mean relative Variance in retention time per peak before alignment:",round(object$summary$variance_before$average,4),"\n"))
    cat(paste("\t","Mean relative Variance in retention time per peak after alignment:",round(object$summary$variance_aligned$average,4),"\n"))

    }


