#' Normalisation of peak abundancies
#'
#' @description
#' Calculates the relative abundance of a peak by normalising an intensity measure with regard to the cumulative abundance of all peaks that are present within an individual sample. The desired measure of peak abundance needs to be included in a column of the input dataset aligned with \code{\link{align_chromatograms}}.
#'
#' @param data
#' Object of class GCalign created with \code{\link{align_chromatograms}} or a list of data frames that contain peak list of individual samples.
#'
#' @param conc_col_name
#' Character giving the name of a column in \code{data} containing a variable describing the abundance of peaks (e.g. peak area or peak height).
#'
#' @param out
#' character defining the format of the returned data. Either "List" or "data.frame".
#'
#' @inheritParams align_chromatograms
#'
#' @return
#' Depending on \code{out} either a list of data frame or a single data frame were rows represent samples and columns relative peak abundances. Abundances are given as percentages.
#'
#'  @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#' @examples
#' ## aligned gc-dataset
#' data("aligned_peak_data")
#' ## returns normalised peak area
#' norm_peaks(data = aligned_peak_data, conc_col_name = "area", rt_col_name = "time")
#'
#' @export
#'
norm_peaks <- function(data, conc_col_name = NULL, rt_col_name = NULL, out = c("data.frame","list")) {

out <- match.arg(out)

## some checks
if (is.null(conc_col_name)) {stop("List containing peak concentration is not specified. Define conc_col_name")}

if (inherits(data, "GCalign")) {
    which <-  "aligned"
    conc_list <- data[[which]][[conc_col_name]]
} else if (is.list(data)) {
    conc_list <- lapply(data, function(fx) fx[[conc_col_name]])
}


# Function to do the calculations
rel_abund <- function(conc_df){
    total_con <- sum(conc_df)
    conc_df <- (conc_df/total_con)
    ## if (percent == TRUE)
    conc_df <- conc_df * 100
    return(conc_df)
}

# Lapply over all data.frames
rel_con_list <- lapply(conc_list, rel_abund)
if (out == "data.frame") {
    if (is.null(rt_col_name)) stop("Specify rt_col_name")
        x <- rel_con_list[-1]
        x <- as.data.frame(t(do.call(cbind,x)))
        colnames(x) <- data[[which]][[rt_col_name]]["mean_RT"][[1]]
        rel_con_list <- x
}
return(rel_con_list)
}
