# Peak Alignment Algorithm to account for shifts in retention of individual
# substances among chromatograms created by GC/GC-MS
# Concept based on Matlab Workflow written by Martin 

# Basic steps
# 1.) loading data from large csv file containing gc data for all individuals from Xcalibur (7 cols per ind)
# 2.) Linear transformation of retention times relative to one reference chromatogram to maximize number of shared peaks
# 3.) Allocation of substances in same rows in all chromatograms based on similar RT?s within a specified interval
# 4.) Repetitive Allocation by randomization of the sample order
# 5.) Merging rows of close (i.e.) duplicated peaks || Not implemented yet!
# 6.) Deleting Substances present in Controls || Not implemented yet!
rm(list=ls())
library(readr)
library(stringr)

source('RetentionCutoff.R')
source('LinearTransformation.R')
source('MatrixOperations.R')
source('CorrectRows.R')
source('ReNaming.R')
source('ChromaVariation.R')
# preprocessing -------------------------------------------------------------------------------

# load data

# gc data
chroma <- read.csv("data/Preen.csv", skip = 4,sep=";") # csv2 introduced mistakes 4 -> 5
chroma <- as.data.frame(apply(chroma, 2, as.numeric)) # convert everything to numeric

# extract names of individuals 
ind_names <- unlist(read_csv2("data/Preen.csv", skip = 1, n_max = 1, col_names = FALSE ))
ind_names <- ind_names[!is.na(ind_names)]
ind_names <- str_replace(ind_names, "-1.raw", "")
ind_names <- tolower(ind_names)

# transform GC data to list of individual GC matrices 
extract <- seq(from = 1, to = ncol(chroma), by = 7)
Chromatograms <- lapply(extract, function(x) chroma[, x:(x+6)])

# name list elements according to individuals
names(Chromatograms) <- ind_names

# rename columns in data.frames
Chromatograms <- lapply(Chromatograms, ReName)


rm(list=c('ind_names','extract','ReName','delete_space_colnames'))
# Start of processing --------------------------------------------------------------------------

# 1.) cut retention times below 8
Chromatograms <- lapply(Chromatograms, RetentionCutoff)

# 2.) Linear Transformation of Retentiontimes
Chroma_aligned <- AlignPeaks(Chromatograms = Chromatograms, References = "w4",
                             Shift=0.05,StepSize=0.005,Error=0.02)

names(Chroma_aligned) <- names(Chromatograms) 

Chromatograms <- Chroma_aligned
# Make List equal in length
Chromatograms <- lapply(Chromatograms, MatrixAppend, Chromatograms)

Length <- (max(sapply(Chromatograms,function(x) out <- MatrixLength(Chromatograms[[1]]))))# To obtain Rows after run of the algorithm
Variation <- mean(VarRows(Chromatograms),na.rm = T)

rm(list=c("Chroma_aligned","chroma","AdjustRetentionTime","AlignPeaks","BestShift",
          "RetentionCutoff","SharedPeaks","PeakShift"))

for (R in 1:1){
  ShuffleOrder <- sample(1:length(Chromatograms))
  Chromatograms <- Chromatograms[ShuffleOrder] # Shuffle
  # Give initial values to some variables
  Error <- 0.021 # allowed deviation around the mean of each row
  Row <- 1
  
  
  while (Row != 'Stop'){ 
    # Loop through all rows, starting with one sample and adding others subsequentially
    # Start with second sample
    # Compare RT with mean of the others
    # Shift current sample down, if RT is above mean
    # Shift all previous samples, if RT is below
    # Leave Samples with RT within error range at their positions
    
    for (S in 2:length(Chromatograms)){                                               
      RT <- Chromatograms[[S]]$ApexRT[Row]                                 
      AvRT <- MeanOfSamples(Chromatograms,Objects =1:(S-1),Row)          
      
      if(is.nan(AvRT)){
        # Do not shift
        AvRT <- RT
      }
      if(RT==0){
        # Do not shift, if current has no substance in that row
        RT <- AvRT
      }
      
      if (RT > AvRT+Error) {                                      
        Chromatograms <- ShiftRows(Chromatograms,S,Row,Error) 
        
      }else if(RT < (AvRT - Error)){
        for(J in 1:(S-1)){     
          # in case that one RT before the current RT is smaller current RT minus Error
          # we think that a second run makes this step unnecessary
          #                 if(Chromatograms[[J]]$ApexRT[Row] < (RT-Error) {
          #                     all_but <- c(1:(S-1))[-J]
          #                     for (i in 1:all_but){
          #                         Chromatograms <- ShiftRows(Chromatograms, i, Row, Error)
          #                     }
          #                 }
          if(Chromatograms[[J]]$ApexRT[Row] <= (RT + Error)){
            # Do nothing, substance?position is valid
          }else{
            Chromatograms <- ShiftRows(Chromatograms,J,Row,Error)
          }
        }
      }
      
    }
    Row <- Row+1
    Chromatograms <- lapply(Chromatograms, MatrixAppend, Chromatograms)# Make all equal in length
    LastSubstance <- max(which(MeanRows(Chromatograms)==max(MeanRows(Chromatograms),na.rm=T)))
    # Remove tail of all-zero rows
    Chromatograms <- lapply(Chromatograms,DeleteZeros) # Remove appended zeros
    
    if (Row>dim(Chromatograms[[S]])[1]){
      Row <- 'Stop' # Signal the end of the iteration
    }
    
  }
  
  Length[R+1]<-max(sapply(Chromatograms,function(x) out <- MatrixLength(Chromatograms[[1]])))
  Variation[R+1] <- mean(VarRows(Chromatograms),na.rm = T)
}

rm(list=c('AvRT','Error','J','Length','R','Row','RT','S','ShuffleOrder','LastSubstance','Variation'))

save(Chromatograms,file="Aligned_Reference_W4_1_Run.RData")


# Post Processing after running the algorithm repeatedly 
# ------------------------------------------------------------------------------------------------------------
rm(list=ls())
source('MatrixOperations.R')
source('MergeRows.R')
source('ChromaVariation.R')
load('Aligned_Reference_W4_1_Run.RData')

# Find rows containing no substances --> Mean of rows == NA, remove them
AverageRTs <- MeanRows(Chromatograms)
if (length(AverageRTs[AverageRTs=="NA"])>0){ #Check for empty rows, and correct accordingly
  Chromatograms <- lapply(Chromatograms,DeleteEmptyRows,AverageRTs)
}

Merging <- 'Start'
while(Merging!='Stop'){
    
  # Find and delete Rows that are 
  # (i) similar in their Retention Time
  # (ii) redundant i.e. no chromatogram contains substances in both rows
  # Always start to pick the fist potential candidate of row-pairs to merge
  # in case merging is not aplicable, take the next position
  # after merging of one pair, update criteria
  # if non is ready to merge anymore, stop the merging algorithm
  
  # Updating merging criterions
  AverageRTs <- MeanRows(Chromatograms) # Average RTs, after merging rows
  similar <- SimilarRows(AverageRTs)    # remaining similarities
  counter <- 1
  while (counter!='Stop'){
    # loop through vector of similar rows, until one shift was done
    Redundant <- sapply(lapply(Chromatograms, CheckRedundancy,similar[counter]),as.vector) #Check first position
    Criterion <- IsRedundant(similar = similar[counter],Redundant=Redundant)
    if(Criterion==1){ # only merge if criterion proves redundancy of one of the rows
      Chromatograms <- lapply(Chromatograms,Merge, ToMerge=similar[counter])
      counter <- 'Stop'
    } else{
      counter <- counter+1
      if (counter>length(similar)){
        Merging <- 'Stop'
      }  
    }
  }
}

save(Chromatograms,file = 'Aligned_Reference_W4_1_Run_Merged_Zeros.RData')

# ----------------------------------------------------------------------------------------------------------------
# Merge close rows, when only a small proportion of samples contain two peaks,
# and prevent the merging in the previous run

Merging <- 'Start'
while(Merging!='Stop'){
  # Find and delete Rows that are 
  # (i) similar in their Retention Time
  # (ii) redundant i.e. no chromatogram contains substances in both rows
  # Always start to pick the fist potential candidate of row-pairs to merge
  # in case merging is not aplicable, take the next position
  # after merging of one pair, update criteria
  # if non is ready to merge anymore, stop the merging algorithm
  
  # Updating merging criterions
  AverageRTs <- MeanRows(Chromatograms) # Average RTs, after merging rows
  similar <- SimilarRows(AverageRTs)    # remaining similarities
  counter <- 1
  while (counter!='Stop'){
    # loop through vector of similar rows, until one shift was done
    Redundant <- sapply(lapply(Chromatograms, CheckRedundancy,similar[counter]),as.vector) #Check first position
    Criterion <- IsRedundant(similar = similar[counter],Redundant=Redundant,Criterion = "Proportional")
    if(Criterion==1){ # only merge if criterion proves redundancy of one of the rows
      Chromatograms <- lapply(Chromatograms,Merge, ToMerge=similar[counter],Criterion="Area")
      counter <- 'Stop'
    } else{
      counter <- counter+1
      if (counter>length(similar)){
        Merging <- 'Stop'
      }  
    }
  }
}

save(Chromatograms,file = 'Aligned_Reference_W4_1_Run_Merged_Proportional.RData')

# ---------------------------------------------------------------------------------------------------------
# Define the Retention times of each chromatogram by the mean retention time among all 
rm(list=ls())
load('Aligned_Reference_W4_1_Run_Merged_Proportional.RData')
source('ChromaVariation.R')
source('CleanChromas.R')
source('MatrixOperations.R')
source('RelativeAbundance.R')
AverageRTs <- MeanRows(Chromatograms)
Chromatograms <- lapply(Chromatograms,EqualRTs,AverageRTs) # Change RTs to mean of given row
Chromatograms <- lapply(Chromatograms,DeleteEmptyRows,AverageRTs) # Remove empty rows
Blanks <- BlankPeaks(Chromatograms)
Chromatograms <- lapply(Chromatograms,DeleteSubstances,Blanks)

SingleChroma <- SingleSubstance(Chromatograms) # Substances present in just one chromatogram
Chromatograms <- lapply(Chromatograms,DeleteSubstances,SingleChroma)

j <- lapply(Chromatograms,ExtractRT)


save(Chromatograms,file = 'Aligned_Cleaned.RData')
#load('Aligned_Cleaned.RData')

Chromatograms <- lapply(Chromatograms,RelativeAbundance)
save(Chromatograms,file='Preen_Scent.RData')
