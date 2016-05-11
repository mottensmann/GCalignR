#'  RT_extract extracts the retention times of all objects of a list
#'  of chromatograms and writes them to a matrix, where rows are samples
#'  and columns retention times. Columnames represent the order of substances
#'  starting with the lowest retention time
#'
#' @param chromatograms List of Chromatograms
#'
#' @return
#'
#'
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @import
#'
#'
#' @export
#'
#'
RT_extract = function(chromatograms,rt_name){

    chromatograms <- lapply(chromatograms,matrix_append,chromatograms) # equal length
    rts <- lapply(chromatograms, function(x) out <- x[1][,rt_name])  # just extract RTs
    rts_matrix <- matrix(NA, nrow = length(rts),ncol = length(rts[[1]]))

    for (C in 1:length(rts)){
        rts_matrix[C,] <- rts[[C]]
    }

    rts_matrix[which(rts_matrix==0)] <- NA # do not consider empty rows

    substances <- colMeans(rts_matrix, na.rm = T)
    ind_names <- names(chromatograms)
    colnames(rts_matrix) <- seq(from=1,to=ncol(rts_matrix))
    rownames(rts_matrix) <- ind_names
    out <- rts_matrix
}
