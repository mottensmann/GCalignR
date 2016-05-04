#' Adds zeros to matrices to make row numbers equal
#'
#' @param Matrix A matrix to append
#' @param Chromatograms List of all matrices
#'
#' @export


matrix_append <- function(Matrix, Chromatograms){
    # Add zeros matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(Chromatograms,function(x) nrow(x)))
    ToAppend <- MaxLength-dim(Matrix)[1]
    Cols <- dim(Matrix)[2]
    Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
    colnames(Zeros) <- names(Matrix)
    Matrix<- rbind(Matrix[,],Zeros)
}


