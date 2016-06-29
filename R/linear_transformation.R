#' Shift peaks to eliminate systematic inaccuracies of peak detection by GC.
#'
#'@description
#'\code{linear_transformation()} applies small linear shifts of all peaks of
#'individual samples with respect to one reference. Thereby the number of shared
#'compounds among samples is maximized. The interval in which linear transformations are evaluated
#'is adjustable as well as the step size within this range.
#'
#'
#' @param reference
#' a character string indicating a sample included in \code{gc_peak_list} used as a reference to align to.
#'
#' @param step_size
#' indicates the step size in which linear shifts are evaluated
#' between \code{max_linear_shift} and \code{-max_linear_shift}.
#'
#' @param error
#' numeric value defining the allowed difference in retention times in
#' derterming if two peaks are shared. The default \code{error=0} counts
#' a peak a shared if retention times match excatly.
#'
#' @inheritParams align_chromatograms
#' @inheritParams matrix_append
#'
#' @return
#' A list of data.frames containing chromatograms with applied linear shifts
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'

linear_transformation <- function(gc_peak_list,reference,
    max_linear_shift=0.05, step_size=0.005, error=0, rt_col_name){
    # This is the master function which calls all sub-functions in order to
    # utilize a maximisation of the number of shared peaks
    # Mandatory arguments of this function are:
    # gc_peak_list = List of Chromatograms, whereby each element of the List is a Matrix with the
    # peak extraction output (7 columns) of Xcalibur
    # References = Name(s) of Reference(s).

    # Include a vector of column names "ColNames" or specifiy the column which holds the
    # Apex of Retention Times "ColumnRT"

    if(file.exists(paste0(as.character(match.call(definition = sys.function(sys.parent(1)), call = sys.call(sys.parent(1)))["data"]),"_LogFile.txt"))){
        sink(paste0(strsplit(as.character(match.call(definition = sys.function(sys.parent(1)), call = sys.call(sys.parent(1)))["data"]),split=".txt"),"_LogFile.txt"),append = TRUE)
        cat("\nSamples in order of comparisons with the reference:\n")
        print(names(gc_peak_list))
        sink()
    }

    ref <- gc_peak_list[[reference]]
    # Chroma_aligned <- list()

    shift_rts <- function(gc_peak_df, ref_df, max_linear_shift, step_size, error) {
        #cat("\n")
        optimal_shift <- peak_shift(gc_peak_df, ref_df, max_linear_shift, step_size, error, rt_col_name)
        shifted <- adjust_retention_time(gc_peak_df, optimal_shift, rt_col_name)
    }

    chroma_aligned <- lapply(gc_peak_list, shift_rts, ref_df = ref, max_linear_shift = max_linear_shift, step_size = step_size, error = error)

    return(chroma_aligned)
}
