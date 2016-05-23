#' Normalisation of peaks
#'
#' @description
#' \code{norm_peaks} calculates the relative abundance of a peaks by normalising
#' with regard to the cumulative abundance of all peaks that are present within an
#' individual chromatogram. The desired measure of peak abundance needs to be a column
#' within the original gas-chromatography data \code{datafile} submitted to
#' \link\code{align_chromatograms}.
#'
#' @param GCout
#' object of class GCaling created with \link\code{align_chromatograms}. Contains a list
#' of data.frames including the retention time and other variables, of which one needs
#' to be named as specified by \code{conc_var_name}.
#'
#' @param conc_var_name
#' character string denoting a column in data.frames of \code{GCout}
#' containing a variable describing the abundance of peaks (e.g. peak area or peak height).
#' @return
#' a list of data.frames containing normalised peak abundances
#' @export
#'
#'

norm_peaks <- function(GCout,conc_var_name=NULL){

    if(class(GCout)!="GCalign"){
        warning("Input is not a output of align_chromatograms, assure the format is correct")
    }
    if(is.null(conc_var_name)){stop("List containing peak concentration is not specified. Define conc_var_name")}
    #########################################################
    # extract concentration measure & convert to a data.frame
    #########################################################

    conc_list <- GCout[["chroma_aligned"]][[conc_var_name]]

    #################################
    # Function to do the calculations
    #################################

    rel_abund <- function(conc_df){
        total_con <- sum(conc_df)
        conc_df <- (conc_df/total_con)*100
        if(sum(conc_df)!=100){
            print('error')
        }
        conc_df
    }
    #############################
    # Lapply over all data.frames
    #############################
    rel_con_list <- lapply(conc_list, rel_abund)
    return(rel_con_list)
}



