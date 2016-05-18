#' Append zeros to matrices
#'
#' @description
#' \code{matrix_append()} adds zeros to a matrix or a data.frames grouped within a list.
#' If required rows containing just zeros are appended to match the dimensions of the
#' largest matrix having the most rows of the list \code{chromatograms}.
#'
#' @details If the matrix has the same dimensions as the largest of \code{chromatograms}
#'          the input matrix is returned.
#'
#'
#' @param matrix a matrix to append.
#'
#' @param chromatograms list of matrices.
#'
#' @return a appended matrix.
#'
#' @export


matrix_append <- function(matrix, chromatograms){
    # Add zeros matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(chromatograms,function(x) nrow(x)))
    ToAppend <- MaxLength-dim(matrix)[1]
    Cols <- ncol(matrix)
    Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
    colnames(Zeros) <- names(matrix)
    matrix<- rbind(matrix[,],Zeros)
}


