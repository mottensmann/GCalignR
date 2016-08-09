#' Summarizing Peak Alignments with GCalignR
#'
#' @description
#' summary method for class "GCalign"
#'
#' @param object
#' an object of class \code{GCalign}, a result of a call to \code{\link{align_chromatograms}}
#'
#' @param write_txt_file
#' logical. If \code{TRUE} the output is passed to \code{.txt} file. No text is printed
#' to the console in this case.
#'
#' @examples
#' summary(gc_peaks_aligned) # prints summary to the Console
#'
#' @param ...
#' optional arguments affecting the summary produced.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @export
#'
summary.GCalign <- function(object,write_txt_file=FALSE,...){

    x <- as.list(match.call())[["object"]] # Name of the GCalign Object

if(write_txt_file==TRUE){
    sink(paste0(x,"_summary.txt"),append = FALSE) # Prepare a Text file
}
    cat("Summary of Peak Alignment running align_chromatograms from package GCalignR\nInput Data: ",
        object[["Logfile"]][["Input"]]["File"])
    cat("\tStart: ", object[["Logfile"]][["Date"]]["Start"],"\tEnd: ",object[["Logfile"]][["Date"]]["End"],"\n\n")
    cat("Number of identified substances among all substances:\n")
    print(object[["Logfile"]][["Aligned"]])
    cat("\n")
    cat(object[["Logfile"]][["Input"]]["Samples"],"Samples were aligned to the reference",paste0("'",object[["Logfile"]][["Input"]]["Reference"],"'"))
    cat(" by means of small linear adjustments\n in order to control for systematic temporal shifts in the gas-chromatography run:\n\n")
    cat("Overview of applied shifts\n")
    print(object[["Logfile"]][["LinearShift"]])
    cat("\nVariation of Retention Times. Minimum and Maximum values refer to the smallest and the largest variation\n")
    cat("within individual peaks respectively.")
    cat("Median, Mean & Standard Deviation summarise the population of peaks.\nNote that a reduction in the variation is expected only for the aligned data!\n\n")
    print(object[["Logfile"]][["Variation"]])

if(write_txt_file==TRUE) {
    sink() # Close the connection
    cat("Summary is saved as",paste0(x,"_summary.txt"))
}
}




