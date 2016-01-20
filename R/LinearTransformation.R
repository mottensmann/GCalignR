AlignPeaks <- function(Chromatograms,References,
    Shift=0.05,StepSize=0.01,Error=0){
    # This is the master function which calls all sub-functions in order to 
    # utilize a maximisation of the number of shared peaks
    # Mandatory arguments of this function are:
    # Chromatograms = List of Chromatograms, whereby each element of the List is a Matrix with the
    # peak extraction output (7 columns) of Xcalibur
    # References = Name(s) of Reference(s).  
    
    # Include a vector of column names "ColNames" or specifiy the column which holds the 
    # Apex of Retention Times "ColumnRT"
 
    Refs <- Chromatograms[References]
    Chroma_aligned <- list()
    for (N in 1:length(Chromatograms)){   
        for (M in 1:length(Refs)){
            Data <- Chromatograms[[N]]
            Ref <- Refs[[M]] 
            OptimalShift <- PeakShift(Data,Ref,Shift,StepSize,Error) #Find the best shift
            Shifted <- AdjustRetentionTime(Data,OptimalShift) # Apply best shift
            
        }
    Chroma_aligned[[N]] <- Shifted
    }
    Chroma_aligned
}

PeakShift <- function(Data,Reference,Shift=0.05,StepSize=0.01,Error=0){
    # This functions shifts retention times of a chromatogram and estimates
    # the number of shared peaks with the reference.
    # Calling 'SharedPeaks' to count the number of shared peaks for a given shift
    # Calling 'BestShift' to find the best shift leading to the maximum similarity
    # If two Shifts lead to the maximu, take the adjustment with the smallest absolut value
    RightShift <- Shift
    LeftShift <- Shift*-1
    ShiftSteps <- seq(from=LeftShift,to=RightShift,by=StepSize)
#     PeaksShared <- rep(0,length(ShiftSteps))
#     PeaksLag <- rep(0,length(ShiftSteps))
    Output <- SharedPeaks(Data,Reference,ShiftSteps,Error) # List containg shared Peaks and their shifts
    Output <- BestShift(Output) # Which is the best setting
    Output # Numeric Value, indicating the best shift (e.g. -0.02 seconds)
}


SharedPeaks <- function(Data,Reference,ShiftSteps,Error=0) {
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

BestShift <- function(List){
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
        if (length(BestFit > 1)) BestFit <- BestFit[1]
    } else{
        BestFit
    }
    
}
AdjustRetentionTime <- function(Data,OptimalShift){
    # Apply the estimated shift of the retention time to the 
    # selected chromatogram to maximize the similarity compared to the reference
    Data$ApexRT <- Data$ApexRT + OptimalShift
    Data
}
