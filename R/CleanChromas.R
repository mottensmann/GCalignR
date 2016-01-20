EqualRTs = function(Chromatograms,AverageRTs){
  # Set Mean Retention Times in every row among chromatograms
Chromatograms$ApexRT <- AverageRTs  
Chromatograms
}

BlankPeaks = function(Chromatograms,Blanks=c("w17","w37","w47","w57","w67","w77")){
  #Select all peak present in Blanks, and delete them in every chromatogram
  Peaks <- numeric(0)
  for (i in 1:length(Blanks)){
    Index <- which(names(Chromatograms)==Blanks[i])
  temp <- as.vector(subset(Chromatograms[[Index]],ApexRT>0 & Area>0)[,'ApexRT'])
  Peaks <- c(Peaks,temp)
  }
  RT <- unique(Peaks)
  Rows <- sapply(RT, function(x) out <- which(Chromatograms[[1]]$ApexRT==x))
  Rows
}

ExtractRT = function(Chromatograms){
  Chromatograms <- Chromatograms$ApexRT
}

DeleteEmptyRows = function(Chromatograms,AverageRTs){
  Rows <- which(!is.na(AverageRTs))
  Chromatograms <- Chromatograms[Rows,]
}


SingleSubstance = function(Chromatograms){
  # Find Substances, which exist only in one chromatogram
  # Therefore first combine everything to big matrix, to select substance column-wise
  Scent <- matrix(NA, nrow = length(Chromatograms),ncol = dim(Chromatograms[[1]])[1])
  for (C in 1:length(Chromatograms)){
    for (R in 1:dim(Chromatograms[[1]])[1]){
      Scent[C,R] <- Chromatograms[[C]]$AreaRel[R]
    }
  }
  
ToDelete <- numeric(0)
  
  for (i in 1:dim(Scent)[2]){
    Substances <- Scent[,i]
    if(length(Substances[Substances>0])<=1){
      ToDelete <- c(ToDelete,i)
    }
  }
    ToDelete
  }

DeleteSubstances = function(Chromatogram, ToDelete){
  Rows <- 1:dim(Chromatogram)[1]
  Valid <- Rows[!Rows %in% ToDelete]
Chromatogram <- Chromatogram[Valid,]
}
  