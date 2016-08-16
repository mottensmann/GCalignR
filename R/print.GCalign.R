#' Summarizing Peak Alignments with GCalignR
#'
#' @description
#' print method for class "GCalign"
#'
#' @param object
#' an object of class \strong{"GCalign"}, the result of a call to \code{\link{align_chromatograms}}
#'
#' @param write_text_file
#' logical. If \code{TRUE} the output is passed to \code{.txt} file. No text is printed
#' to the console in this case.
#'
#' @examples
#' print(gc_peaks_aligned) # prints summary to the Console
#'
#' @param ...
#' optional arguments passed on to methods. Internally \link{cat} is used to generate output, hence optional arguments are currently not ignored.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @export
#'
print.GCalign <- function(x,write_text_file=FALSE,...){

    # Extract some informative pieces from the Logfile
    object_name <- as.list(match.call())[["object"]] # Name of the GCalign Object
    function_call <- x[["Logfile"]][["Call"]][1,]
    function_call[["sep"]] <- ifelse(function_call[["sep"]]=="\t","\\t",function_call[["sep"]])
    names_call <- names(function_call)
    peaks <- x[["Logfile"]][["Aligned"]]
    names_peaks <- names(peaks)

if(write_text_file==TRUE){
    sink(paste0(object_name,"_summary.txt"),append = FALSE) # Prepare a Text file
}
    cat(stringr::str_wrap(paste("Summary of Peak Alignment running align_chromatograms from package GCalignR\nInput Data: ",x[["Logfile"]][["Input"]]["File"]),width=80,exdent=2,indent=2))
    cat("\tStart: ", x[["Logfile"]][["Date"]]["Start"],"\tEnd: ",x[["Logfile"]][["Date"]]["End"],"\n\n")

    cat("Call:\n")
    cat("GCalignR::align_chromatograms(")

    c <- 1
    text <- character()
    for(i in names_call){
        if(c < length(names_call)){
     text <- c(text,paste0(i,"=",as.character(function_call[[i]]),""))
            }else{
     text <- c(text,paste0(i,"=",as.character(function_call[[i]]),")"))
        }
        c <- c +1
    }
    cat(stringr::str_wrap(paste(text,collapse = ", "),width=80,exdent=2,indent=2))


    cat("\n\nSummary of scored substances:\n\n")
    cat("In total",peaks["Peaks"],"substances were identified among all samples.\n")
    if(any((names_peaks) %in% "In_Blanks")) cat(peaks["In_Blanks"],"substances were present in blanks. The corresponding peaks as well as the blanks were removed from the data set.")
    if(any((names_peaks) %in% "Singular")) cat(peaks["Singular"],"substances were present in just one single sample and were removed.")
    cat("After the processing",peaks["Retained"], "substances are reatained in the data set")
    cat("\n")
    cat(x[["Logfile"]][["Input"]]["Samples"],"Samples were aligned to the reference",paste0("'",x[["Logfile"]][["Input"]]["Reference"],"'"))
    cat(" by means of small linear adjustments\n in order to control for systematic temporal shifts in the gas-chromatography run:\n\n")
   # cat("Overview of applied shifts\n")
    #print(x[["Logfile"]][["LinearShift"]])
    cat("\nVariation of Retention Times. Minimum and Maximum values refer to the smallest and the largest variation\n")
    cat("within individual peaks respectively.")
    cat("Median, Mean & Standard Deviation summarise the population of peaks.\nNote that a reduction in the variation is expected only for the aligned data!\n\n")
    print(x[["Logfile"]][["Variation"]])

if(write_text_file==TRUE) {
    sink() # Close the connection
    cat("Summary is saved as",paste0(object_name,"_summary.txt"))
}
}




