#' print.GCalign: Prints a the protocol of an alignment process exectuted via align_chromatogtams
#'
#' @details
#' Prints the LogFile written during the alignment-process, if one was created.
#'
#' @param object
#' Object of class \code{GCalign}
#' @param ... additional input to summary.GCalign
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'  @export
summary.GCalign <- function(object,...){

    if(file.exists(paste0(strsplit(object$call$data,split = ".txt"),"_LogFile.txt"))){
        text<-readLines(paste0(strsplit(object$call$data,split = ".txt"),"_LogFile.txt"))
    cat(text,sep = "\n")
    }else{
        cat("No LogFile was created!","\nSet LogFile=TRUE when executing align_chromatograms()")
    }
}




