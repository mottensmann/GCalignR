#' Normalisation of peaks
#'
#' @description
#' \code{norm_peaks} calculates the relative abundance of a peak by normalising
#' with regard to the cumulative abundance of all peaks that are present within an
#' individual chromatogram. The desired measure of peak abundance needs to be a column
#' within the original gas-chromatography dataset aligned by \link{align_chromatograms}.
#'
#' @param data
#' Object of class GCalign created with \link{align_chromatograms}. Contains a list
#' of data frames including the retention time and other variables, of which one needs
#' to be specified as \code{conc_col_name}.
#'
#' @param conc_col_name
#' Character string denoting a column in data frames of \code{data}
#' containing a variable describing the abundance of peaks (e.g. peak area or peak height).
#'
#' @param out
#' character string defining the format of the returned data. Either "List" or "data.frame".
#'
#' @inheritParams align_chromatograms
#'
#' @return
#' Depending on \code{out} either a list of data frame or a single data frame were rows represent samples and columns relative peak abundancies. Abundancies are given in percent.
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
which <- "aligned" # which <- match.arg(which)

## some checks
if (class(data) != "GCalign") {warning("Input is not a output of align_chromatograms, assure the format is correct")}
if (is.null(conc_col_name)) {stop("List containing peak concentration is not specified. Define conc_col_name")}

conc_list <- data[[which]][[conc_col_name]]

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
## unused code chunks
## ----------------------------------------------
# if (which == "raw") which <- "input_matrix"
## additional function parameters
#  ..., percent = TRUE, which = c("aligned","raw")
# param percent
# By default percent values are returned (i.e. relative abundancies scale up to 100.) If
#
# param which
# character string naming the data source to normalise. Either the aligned data or the raw data prior to alignment. The latter can be used for demonstration purposes.
# -----------------------------------------------
}
