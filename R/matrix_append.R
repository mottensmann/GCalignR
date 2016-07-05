#' Append zeros to matrices
#'
#' @description
#' \code{matrix_append()} adds zeros to a a data.frames grouped within a list.
#' If required rows containing just zeros are appended to match the dimensions of the
#' largest matrix having the most rows of the list \code{gc_peak_list}.
#'
#' @details If the matrix has the same dimensions as the largest of \code{gc_peak_list}
#'          the input matrix is returned.
#'
#'
#' @param gc_peak_df
#' data.frame containing GC-data (e.g. retention time, peak area, peak height) of one sample. Variables are stored in columns.
#' Rows represent distinct peaks. Retention time is a required variable.
#'
#' @param gc_peak_list
#' List of data.frames. Each data.frame contains GC-data (e.g. retention time, peak area, peak height) of one sample. Variables are stored in columns.
#' Rows represent distinct peaks. Retention time is a required variable.
#'
#' @return a appended matrix.
#' @keywords internal
#' @export
#'
#'

matrix_append <- function(gc_peak_df, gc_peak_list){
    # Add zeros matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(gc_peak_list,function(x) nrow(x)))
    ToAppend <- MaxLength-nrow(gc_peak_df)
    Cols <- ncol(gc_peak_df)
    Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
    colnames(Zeros) <- names(gc_peak_df)
    gc_peak_df<- rbind(gc_peak_df[,],Zeros)
    return(gc_peak_df)
}


