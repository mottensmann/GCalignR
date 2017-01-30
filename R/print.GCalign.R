#' Summarizing Peak Alignments with GCalignR
#'
#' @description
#' print method for class "GCalign"
#'
#' @param x
#' an object of class \strong{"GCalign"}, the result of a call to \code{\link{align_chromatograms}}
#'
#' @param write_text_file
#' logical. If \code{TRUE} the output is passed to \code{.txt} file. No text is printed
#' to the console in this case.
#'
#' @examples
#' print(aligned_peak_data) # prints summary to the Console
#'
#' @param ...
#' optional arguments passed on to methods. Internally \link{cat} is used to generate output, hence optional arguments are currently not ignored.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#' @import stringr
#' @export
#'
print.GCalign <- function(x,write_text_file=FALSE,...){

    object_name <- as.list(match.call())[["x"]] # Name of the GCalign Object

    # Extract some informative pieces from the Logfile
    function_call <- x[["Logfile"]][["Call"]] # List of all function arguments
    function_call[["sep"]] <- ifelse(function_call[["sep"]]=="\t","\\t",function_call[["sep"]])
    function_call[["blanks"]] <- if(length(function_call[["blanks"]])>2){
        paste0("(",paste(as.character(function_call[["blanks"]][-1]),collapse = ", "),")")
    }
    function_call[["write_output"]] <- if(length(function_call[["write_output"]])>1){
        paste0("(",paste(as.character(function_call[["write_output"]][-1]),collapse = ", "),")")

    }
    # Tweak the NULL entries, so they are printable
    which_NULL <- as.vector(which(lapply(function_call, function(x) out <- class(x))=="NULL"))
    function_call[which_NULL] <- "NULL"
    names_call <- names(function_call)
    peaks <- x[["Logfile"]][["Aligned"]]
    names_peaks <- names(peaks)

if(write_text_file==TRUE){ # write a textfile, no output to the Console if TRUE
    sink(paste0(object_name,"_summary.txt"),append = FALSE) # Prepare a Text file
}

# Information about the date
# --------------------------
    cat(stringr::str_wrap(paste("Summary of Peak Alignment running align_chromatograms from package GCalignR\nInput: ",x[["Logfile"]][["Input"]]["File"]),width=80,exdent=2,indent=2))
    cat("\tStart: ", x[["Logfile"]][["Date"]]["Start"],"\tFinished: ",x[["Logfile"]][["Date"]]["End"],"\n\n")
    cat("Call:\n")

# Function Call
# ---------------
    c <- 1
    text <- character()# used in this function to collect and then print output
    text <- c(text, (paste0("GCalignR::align_chromatograms(")))
    for(i in names_call){
        if(c < length(names_call)){
     text <- c(text,paste0(i,"=",paste0(as.character(function_call[[i]]),", ")))
            }else{
     text <- c(text,paste0(i,"=",as.character(function_call[[i]]),")"))
        }
        c <- c +1
    }
    cat(stringr::str_wrap(paste(text,collapse = ""),width=80,exdent=2,indent=2))

    # Summary of substances that have been scored and filtred
# ---------------------------------------------------------
text <- character()
    cat("\n\nSummary of scored substances:\n\n")
print(peaks)
cat("\n")
text <- c(text,paste("In total",peaks["Peaks"],"substances were identified among all samples."))

if(any((names_peaks) %in% "In_Blanks")) text <- c(text,paste(peaks["Blanks"],"substances were present in blanks. The corresponding peaks as well as the blanks were removed from the data set."))
    if(any((names_peaks) %in% "Singular")) text <- c(text,paste(peaks["Singular"],"substances were present in just one single sample and were removed."))
text <- c(text,paste(peaks["Retained"], "substances are retained after all filtering steps."))
cat(stringr::str_wrap(paste(text,collapse = " "),width=80,exdent=2,indent=2))

# Reference & Sample  Names.
# -------------------
text <- character()
cat("\n\nSample Overview")
text <- c(text,paste("The following",x[["Logfile"]][["Input"]]["Samples"],"Samples were aligned to the reference",paste0("'",x[["Logfile"]][["Input"]]["Reference"],"':")))
cat(stringr::str_wrap(paste(text,collapse = " "),width=80,exdent=2,indent=2))
cat("\n")
text <- paste(names(x[["aligned"]][[1]])[-1],collapse = ", ") # Sample Names
cat(stringr::str_wrap(paste(text,collapse = " "),width=80,exdent=2,indent=2))
cat("\n\n")

# Refer to plots
# -------------
cat("For further details:\n")
text <- paste0("Type 'gc_heatmap(",object_name,")'"," to retrieve a heatmap for the alignment accuracy")
cat(stringr::str_wrap(paste(text,collapse = " "),width=80,exdent=2,indent=2))
cat("\n")
text <- paste0("Type 'plot(",object_name,")'"," to retrieve further diagnostic plots")
cat(stringr::str_wrap(paste(text,collapse = " "),width=80,exdent=2,indent=2))

if(write_text_file==TRUE) {
    sink() # Close the connection
    cat("Summary is saved as",paste0(object_name,"_summary.txt"))
}
}




