#' Extract retention times to plot a heatmap
#'
#' @description \code{rt_extract()} extracts retention times for all samples within a
#'              a list of chromatograms and formats them to one matrix.
#'
#' @inheritParams matrix_append
#'
#' @inheritParams align_chromatograms
#'
#'
#' @return data.frame containing retention times of all samples in \code{gc_peak_list}.
#'          Samples are sorted in rows, columns represent substances. Columns are named
#'          by the mean retention time based on all samples containing a retention time
#'          > 0 at this location.
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'
rt_extract <- function(gc_peak_list,rt_col_name){

    # blanks and del_single_sub are removed, since their removal
    # is of importance only for the last step, where it is applied
    # outside this function call

    ###############################################
    # Make length equal, if differences are present
    ###############################################

    gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)
    id <- names(gc_peak_list)

    ####################################################################
    # optional, depends on arguments regarding blanks and del_single_sub
    ####################################################################

    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))


    #################################
    # calculate final retention times
    #################################

    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))
    rt_mat <- as.data.frame(t(rt_mat))
    rt_mat2 <- rt_mat
    rt_mat2[rt_mat2==0] <- NA
    colnames(rt_mat) <-
        as.character(round(colMeans(rt_mat2,na.rm = T),3))
    rt_mat <- cbind(id,rt_mat)



}
