Merge = function(Chromatograms,ToMerge,Criterion="Zero"){
  # Check always the row containing just zeros, in case of zeros in both, just delete one of them
  # To Merge == Last row of a similar pair
    Row1 <- ToMerge-1 # Previous Row
    Row2 <- ToMerge # Current Row
    R1 <- Chromatograms$ApexRT[Row1]
    R2 <- Chromatograms$ApexRT[Row2]
    if (Criterion=="Zero"){
    if (Row1>1 & Row2<dim(Chromatograms)[1]){ # Avoid taking first and last rows
    if(R1==0){ 
      #  Delete Row1
      Chromatograms <- rbind(Chromatograms[1:(Row1-1),],Chromatograms[Row2:dim(Chromatograms)[1],])
    } else if(R2==0){
      # Delete Row2
      Chromatograms <- rbind(Chromatograms[1:Row1,],Chromatograms[(Row2+1):dim(Chromatograms)[1],])
    }
    }
    }
    if(Criterion=="Area"){ # Take the larger of two peaks, defined by the are of the peaks 
      if (Chromatograms$Area[Row1]>=Chromatograms$Area[Row2]){
        Chromatograms <- rbind(Chromatograms[1:Row1,],Chromatograms[(Row2+1):dim(Chromatograms)[1],])
      } else if (Chromatograms$Area[Row1]<Chromatograms$Area[Row2]) {
        Chromatograms <- rbind(Chromatograms[1:(Row1-1),],Chromatograms[Row2:dim(Chromatograms)[1],])
        
        }
    }
    Chromatograms
    }
    
    


