#' Normalisation of peaks
#'
#' @description
#' \code{norm_peaks} calculates the relative abundance of a peaks by normalising
#' with regard to the cumulative abundance of all peaks that are present within an
#' individual chromatogram. The desired measure of peak abundance needs to be a column
#' within the original gas-chromatography data \code{datafile} submitted to
#' \link{align_chromatograms}.
#'
#' @param GCout
#' object of class GCaling created with \link{align_chromatograms}. Contains a list
#' of data.frames including the retention time and other variables, of which one needs
#' to be named as specified by \code{conc_col_name}.
#'
#' @param conc_col_name
#' character string denoting a column in data.frames of \code{GCout}
#' containing a variable describing the abundance of peaks (e.g. peak area or peak height).
#'
#' @param percent
#' logical, if \code{percent = TRUE} normalised conentraion scale up to 100, otherwise to 1.
#' @return
#' a list of data.frames containing normalised peak abundances
#'
#'  @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @keywords internal
#' @export
#'
norm_peaks <- function(GCout,conc_col_name = NULL, rt_col_name = NULL, out=c("list","data.frame"), percent = TRUE){
out <- match.arg(out)

# some checks
if (class(GCout) != "GCalign") {warning("Input is not a output of align_chromatograms, assure the format is correct")}
if (is.null(conc_col_name)) {stop("List containing peak concentration is not specified. Define conc_col_name")}

conc_list <- GCout[["aligned"]][[conc_col_name]]

# Function to do the calculations
rel_abund <- function(conc_df){
    total_con <- sum(conc_df)
    conc_df <- (conc_df/total_con)
    if (percent == TRUE) conc_df <- conc_df * 100
    return(conc_df)
}

# Lapply over all data.frames
rel_con_list <- lapply(conc_list, rel_abund)
if (out == "data.frame") {
        x <- rel_con_list[-1] #
        x <- as.data.frame(t(do.call(cbind,x)))
        colnames(x) <- GCout[["aligned"]][[rt_col_name]]["mean_RT"][[1]]
        rel_con_list <- x
}
return(rel_con_list)
}
