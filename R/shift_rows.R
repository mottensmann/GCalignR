#' inserts row of zeros in chromatogram
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param current_sample_index index of sample to insert zeros
#' @param retention_row current retention time row in which to insert zeros
#'
#' @return
#' chromatogram \code{data.frame} with inserted zeros in a given row
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @keywords internal
#'@export

shift_rows = function(chromatograms, current_sample_index, retention_row){
    # Shift Rows to sort substances by rows
    # List of Chromatograms
    # Index is the position of the Object in Chromatogram, which needs to be shifted
    n_col<- ncol(chromatograms[[1]])
    zeros <- as.data.frame(matrix(0,nrow=1,ncol=n_col))           # Create zeros
    colnames(zeros) <- names(chromatograms[[1]])

    chroma_temp <-  chromatograms[[current_sample_index]]

    if(retention_row != 1){
        chroma_temp<- rbind(chroma_temp[1:(retention_row-1), ], zeros,
                            chroma_temp[retention_row:nrow(chroma_temp), ])
    } else {
        chroma_temp <- rbind(zeros, chroma_temp)
    }
    # submit the new list, where the shift was applied to the Object at position current_sample_index
    chromatograms[[current_sample_index]] <- chroma_temp
    chromatograms
}
