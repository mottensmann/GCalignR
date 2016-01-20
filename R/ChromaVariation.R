# Evaluate Algorithm performance
MeanRows = function(Chromatograms){
  MeanRow <- numeric(0)
  for (N in 1:max(sapply(Chromatograms,MatrixLength))){ #Loop through all rows
    MeanRow[N]<-MeanOfSamples(Chromatograms,c(1:length(Chromatograms)),N)
  }
  MeanRow
}

VarRows = function(objects){
  SD <- numeric(0)
  for (N in 1:max(sapply(Chromatograms,MatrixLength))){
    SD[N]<-VarOfSamples(Chromatograms,c(1:length(Chromatograms)),N)
  }
  SD
}

VarOfSamples = function(Data,Objects, Row){
  # Estimate the Variation in Retention Times within Rows  
  # NA indicates that only one Substance exists  
  SD <- numeric(0)
  for (N in 1:length(Objects)){
    SD <- c(SD,Data[[Objects[N]]]$ApexRT[Row])
  }  
  SD <- sd(SD[SD>0],na.rm = T)       # Neglect Zeros (i.e. do not incorporate samples without the substance in the Row)
  SD
}


SimilarRows = function(AverageRTs,MinDistance=0.05){
  # Estimate Degree of Similarity between subsequent rows by comparison of mean
  # retention times of rows.
  # Similarity is evaluated at the level of MinDistance, given in seconds
  # The output gives the the position of a row, that is similar to the previous
  # i.e. row 5 is similar to 4
  
  Difference <- rep(NA, (length(AverageRTs)-1)) # Estimate Difference between adjacent rows
                    for (i in 2:length(AverageRTs)){
                      Difference[i] <- AverageRTs[i]-AverageRTs[i-1]  
                    }
                    similar <- which(Difference<=MinDistance) # Which rows differ less the MinDistance ?
                    similar
}
                     
CheckRedundancy = function(Chromatograms,similar){
  # If only one of two neighbouring rows contain a substance
  # they are redundant, coded by a One
  Row1 <- Chromatograms$ApexRT[similar-1] # Extract previous row
  Row2 <- Chromatograms$ApexRT[similar] # Extract current row
  Redundant <- 0
  if (Row1==0 | Row2==0){
    Redundant <- 1
  }
  Redundant
}

IsRedundant = function(similar, Redundant, Criterion="Strict"){
  # Indicates by a binary output variable (1/0) if rows should be merged
  # Methods: Strict: A single sample with two peaks prevents merging
  #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
  ToMerge <- 0
      if (Criterion=="Strict"){
      if(sum(Redundant)/length(Redundant)==1){
      ToMerge <- 1
      }
  } else if (Criterion=="Proportional"){
    if(sum(Redundant)/length(Redundant)>=0.95){
      ToMerge <- 1
    }
  }
  ToMerge
  }

MeanOfSamples = function(Data,Objects, Row){
  # Calculate mean retention times within a row for a number of objects
  # Data is list of Chromatograms
  # Objects is a vector specifying the objects on which the calculation of the mean is based
  # Row indicates the row to look at
  
  AvRT <- numeric(0)
  for (N in 1:length(Objects)){
    AvRT <- c(AvRT,Data[[Objects[N]]]$ApexRT[Row])
  }  
  AvRT <- mean(AvRT[AvRT>0])       # Neglect Zeros (i.e. do not incorporate samples without the substance in the Row)
  AvRT <- round(AvRT,digits = 2)
}