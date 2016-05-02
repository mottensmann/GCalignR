

ShiftRows = function(Chromatograms,Index,Row,Error=0.02){
    # Shift Rows to sort substances by rows
    # List of Chromatograms
    # Index is the position of the Object in Chromatogram, which needs to be shifted
    Cols<- dim(Chromatograms[[1]])[2]
    Zeros <- as.data.frame(matrix(0,nrow=1,ncol=Cols))           # Create zeros
    colnames(Zeros) <- names(Chromatograms[[1]])
    
    Chroma_temp <-  Chromatograms[[Index]] 
    
    if(Row!=1){
        Chroma_temp<- rbind(Chroma_temp[1:(Row-1),],Zeros,Chroma_temp[Row:dim(Chroma_temp)[1],])
    } else {
        Chroma_temp <- rbind(Zeros, Chroma_temp)
    }
    Chromatograms[[Index]] <- Chroma_temp # submit the new list, where the shift was applied to the Object at position Index
    Chromatograms
}  

