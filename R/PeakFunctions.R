DataSorter=function(DataFile, DataNames, DataCols=7, SkipRows=5,NewFolder) {
  # Split the output of Xcalibur into single files for each sample
  # DataFile: CSV-File containing data belonging to many samples, 
  # tables horizontally structured
  # DataNames = CSV-File holding names of samples in a single column
  # NOTE: The order needs to be the same as in the Xcalibur-File
  # DataCols specifies how many columns belong to each sample, By default its 7
  # Skip rows indicates how many rows of describtive text are above the data table
  # Note: Names of columns are dropped
  
  data <- read.csv(DataFile,header=F, skip=SkipRows,sep=";")
  names <- read.csv(DataNames,header=F,sep=";")
  rows <- dim(data)[1]
  

  if (!file.exists(NewFolder)){
    dir.create(file.path(getwd(),NewFolder)) # Create a folder for individual files
  }else{
    warning('Data was already splitted for each individual')
  }
  OldDir <- getwd()
  setwd(NewFolder)
  
  
  for (N in 1:dim(names)[1]){                 # Loop through the Vector of Samples
    first <- (N*DataCols)-DataCols+1          # Select the corresponding columns
    last <- first+DataCols-1
    
    temp <- data[,first:last]
    write.csv(temp,file=paste0(as.character(names[N,]),".csv"))
  }
  setwd(OldDir)
  Output <- names
}

AlignPeaks = function(Chromatograms,References,ColNames=NULL,ColumnRT=NULL,
                      MinRetentionTime=8,Shift=0.05,StepSize=0.01,Error=0){
  # Tis is the master function which calls all sub-functions in order to 
  # utilize a maximisation of the number of shared peaks
  # Mandatory arguments of this function are:
  # Chromatograms = vector containing the names (without .csv) of single semicolon-delimited csv-files
  # References = vector containing the name of one Reference.  
  
  # Include a vector of column names "ColNames" or specifiy the column which holds the 
  # Apex of Retention Times "ColumnRT"
  for (N in 1:length(Chromatograms)){   
    for (M in 1:length(References)){
      Data <- FileImport(Chromatograms[N],ColNames,ColumnRT,Folder="Preen_Samples")   #Read the data
      Reference <- FileImport(References[M],ColNames,ColumnRT)
      
      Data <- RetentionCutoff(Data,Low=MinRetentionTime)           #Cut-off small retention times
      Reference <- RetentionCutoff(Reference,Low=MinRetentionTime)
      
      OptimalShift <- PeakShift(Data,Reference,Shift,StepSize,Error) #Find the best shift
      Shifted <- AdjustRetentionTime(Data,OptimalShift) # Apply best shift
      SaveShiftedChromas(Shifted,Chromatograms,N) # Save the alligned Chromatogram in a Folder called adjusted
    }}
  
}

FileImport = function(ID, ColNames=NULL,ColumnRT=1,Folder="Adjusted"){ # Change to Adjut again !!!
  # Import a Chromatogram by submitting the ID of a sample (e.g. "W32")
  # Loads the corresponding csv-file (e.g. "W32.csv")
  # ColNames or the Column holding the retention times must be specified
  data <- read.csv(paste0(Folder,"/",as.character(ID),".csv"))
  data <- data[,2:8]
  if (!is.null(ColNames)){
    colnames(data) <- ColNames
  } else if (is.null(ColNames) & (is.null(ColumnRT))) {
stop('Column containing the Retention time is not specified !')      
  } else {
    colnames(data)[ColumnRT] <- "ApexRT"
  }
  data <- subset(data,ApexRT != 'NA') # Remove NAs at the end 
      }
RetentionCutoff = function(Data, Low=8, High=NULL){
  # Removes all Retentionstimes that are below and/or above the
  # given Threshold. By Default there is only a lower boundary
if (is.null(High)){
  temp <- subset(Data, ApexRT > Low)
}else{
  temp <- subset(Data, ApexRT > Low & ApexRT< High)
}
  Data <- temp
}

PeakShift = function(Data,Reference,Shift=0.05,StepSize=0.01,Error=0){
  # This functions shifts retention times of a chromatogram and estimates
  # the number of shared peaks with the reference.
  # Calling 'SharedPeaks' to count the number of shared peaks for a given shift
  # Calling 'BestShift' to find the best shift leading to the maximum similarity
  # If two Shifts lead to the maximu, take the adjustment with the smallest absolut value
  RightShift <- Shift
  LeftShift <- Shift*-1
  ShiftSteps <- seq(from=LeftShift,to=RightShift,by=StepSize)
  PeaksShared <- rep(0,length(ShiftSteps))
  PeaksLag <- rep(0,length(ShiftSteps))
  Output <- SharedPeaks(Data,Reference,ShiftSteps,Error) # List containg shared Peaks and their shifts
  Output <- BestShift(Output) # Which is the best setting
  Output # Numeric Value, indicating the best shift (e.g. -0.02 seconds)
}

SharedPeaks = function(Data,Reference,ShiftSteps,Error=0) {
  # Calculate the Number of shared peaks between a Chromatogram and its reference
  # Shared peaks fall within the retention time of the reference and +- the Error [s]
  Ref <- Reference$ApexRT
  NoOfPeaks <- numeric(0)
  for (j in 1:length(ShiftSteps)){
    Temp <- Data$ApexRT+ShiftSteps[j] # Shift all Peaks by the same step 
    Peaks <- 0
    for (k in 1:length(Temp)){ # loop through all Peaks of the current Sample
      for (l in 1:length(Ref)){ # loop through the Reference Chromatogram
        TempPeak <- Temp[k]
        RefPeak <- Ref[l]
        if (TempPeak!=0){ # Avoid comparison with cases of RT=0
        if ((TempPeak <= RefPeak+Error) & (TempPeak>=RefPeak-Error)){
        Peaks=Peaks+1
        }
        }
      }
      
    }
    NoOfPeaks <- c(NoOfPeaks,Peaks)
    
  }
  Output <- list(NoOfPeaks,ShiftSteps)
}

BestShift=function(List){
  # Selects the optimal Shifting time leading to the maximum number of shared peaks
  # If competing shifts exists (i.e. same number of peaks are shared)
  # selects the smallest absolute value of a shift
  ShiftTime <- as.vector(List[[1]])
  Peaks <- as.vector(List[[2]])
  Index <- which(ShiftTime==max(ShiftTime)) # find the best fit = maximum of shared peaks
  BestFit <- Peaks[Index]
  BestFit
  if (length(BestFit)>1){
    temp <- min(abs(BestFit)) # If equal shifts were found, take the smallest shift applied
  BestFit <- BestFit[BestFit==temp | BestFit==temp*-1]
  if (length(BestFit)>1){
    BestFit <- BestFit[1] # Take only one estimate if there are still two 
  }
  } else{
    BestFit
  }
  
}

AdjustRetentionTime = function(Data,OptimalShift){
  # Apply the estimated shift of the retention time to the 
  # selected chromatogram to maximize the similarity compared to the reference
Data$ApexRT <- Data$ApexRT + OptimalShift
Data
}

SaveShiftedChromas = function(Data,Chromatograms,N){
  # Save the shifted Chromatogram in the Folder Adjusted
  # Data= Data frame holding the chromatogram
  # Chromatogramsis a vector of sample IDs
  # N, the current position in Chromatograms
  if (!file.exists(paste0(getwd(),"/Adjusted"))){
    dir.create(file.path(getwd(),"Adjusted"))
    }
  OldDir <- getwd()
  setwd('Adjusted')
  write.csv(Data,file=paste0(as.character(Chromatograms[N]),".csv"))
  setwd(OldDir)
}

RandomizeOrder = function(Samples){
  # Shuffle the order of Chromatograms randomly
  N <- length(Samples) # How many samples are there
  RandOrder <- sample(1:N)
  NewOrder <- character(0)
  for (i in 1:N){
    NewOrder[i] <- as.character(Samples[RandOrder[i]])
  }
  NewOrder
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

ShiftRows = function(Chromatograms,Index,Row,Error=0.02){
  # Shift Rows to sort substances by rows
  # List of Chromatograms
  # Index is the position of the Object in Chromatogram, which needs to be shifted
  Cols<- dim(Chromatograms[[1]])[2]
  Zeros <- as.data.frame(matrix(0,nrow=1,ncol=Cols))           # Create zeros
  colnames(Zeros) <- names(Chromatograms[[1]])
      if(Row!=1){
      Chromatograms[[Index]] <- rbind(Chromatograms[[Index]][1:(Row-1),],Zeros,Chromatograms[[Index]][Row:dim(Chromatograms[[Index]])[1],])
       
    } else {
      Chromatograms[[Index]] <- rbind(Zeros,Chromatograms[[Index]])
    }
  Chromatograms # submit the new list, where the shift was applied to the Object at position Index
  }  
      

MatrixLength = function(Matrix){
  # Count the number of Rows within a matrix
  dim(Matrix)[1]
}

MatrixAppend = function(Matrix){
  # Add zeros matrices to fit the dimensions of the largest matrix 
  MaxLength <- max(sapply(Chromatograms,MatrixLength))
  ToAppend <- MaxLength-dim(Matrix)[1]
  Cols <- dim(Matrix)[2]
  Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
  colnames(Zeros) <- names(Matrix)
Matrix<- rbind(Matrix[,],Zeros)
}      
    