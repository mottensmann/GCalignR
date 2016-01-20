MatrixLength <- function(Matrix){
    # Count the number of Rows within a matrix
    dim(Matrix)[1]
}

MatrixAppend <- function(Matrix, Chromatograms){
    # Add zeros matrices to fit the dimensions of the largest matrix 
    MaxLength <- max(sapply(Chromatograms,MatrixLength))
    ToAppend <- MaxLength-dim(Matrix)[1]
    Cols <- dim(Matrix)[2]
    Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
    colnames(Zeros) <- names(Matrix)
    Matrix<- rbind(Matrix[,],Zeros)
}   

DeleteZeros = function(Chromatograms){
  
  Chromatograms <- subset(Chromatograms[1:LastSubstance,])
  
}

