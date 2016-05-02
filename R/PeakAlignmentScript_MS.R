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

# source('R/RetentionCutoff.R')
# source('R/LinearTransformation.R')
# source('R/MatrixOperations.R')
# source('R/CorrectRows.R')
# source('R/ReNaming.R')
# source('R/ChromaVariation.R')
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
# extract <- seq(from = 1, to = ncol(chroma), by = 7)
# Chromatograms <- lapply(extract, function(x) chroma[, x:(x+6)])

# name list elements according to individuals
# names(Chromatograms) <- ind_names

# rename columns in data.frames
# Chromatograms <- lapply(Chromatograms, ReName)


# rm(list=c('ind_names','extract','ReName','delete_space_colnames'))

# extract names of individuals 
# ind_names <- unlist(read_csv2("data/Preen.csv", skip = 1, n_max = 1, col_names = FALSE ))
# ind_names <- ind_names[!is.na(ind_names)]
# ind_names <- str_replace(ind_names, "-1.raw", "")
# ind_names <- tolower(ind_names)

# start of package?

# transform GC data to list of individual GC matrices 
# GC_mat_to_list <- function(all_gc_mat, ind_names, var_names) {
#     extract <- seq(from = 1, to = ncol(all_gc_mat), by = length(var_names))
#     chromatograms <- lapply(extract, function(x) chromatograms[, x:(x+length(var_names))])
#     names(chromatograms) <- ind_names
#     chromatograms <- lapply(chromatograms, ReName)
# }


source("R/conv_gc_mat_to_list.R")
source("R/rename_cols.R")
chromatograms <- conv_gc_mat_to_list(chroma, ind_names, var_names = c("RT", "startRT", "endRT", "Area", "xArea", "Height", "xHeight"))


chromatograms <- chromatograms[1:10]
# rm(list=c('ind_names','extract','ReName','delete_space_colnames'))
# Start of processing --------------------------------------------------------------------------

source("R/rt_cutoff.R")
# 1.) cut retention times below 8
chromatograms <- lapply(chromatograms, rt_cutoff, low = 11, high = 20, rt_col_name = "RT")

source("R/linear_transformation.R")
# 2.) Linear Transformation of Retentiontimes

## thinking about reference: default is chromatogram with most peaks - optional: manual 
chroma_aligned <- linear_transformation(chromatograms, shift=0.05, step_size=0.01, error=0, reference = "w3", rt_col_name = "RT")

# Make List equal in length
source("R/matrix_append.R")
chromatograms <- lapply(chroma_aligned, matrix_append, chroma_aligned)


source("R/evaluate_chroma.R")
Length <- (max(unlist(lapply(chromatograms, function(x) out <- nrow(x))))) # To obtain Rows after run of the algorithm
Variation <- mean(var_per_row(chromatograms),na.rm = T)

#rm(list=c("Chroma_aligned","chroma","AdjustRetentionTime","AlignPeaks","BestShift",
#          "RetentionCutoff","SharedPeaks","PeakShift"))

source("R/align_individual_peaks.R")
source("R/shift_rows.R")
chromatograms_aligned <- align_individual_peaks(chromatograms, error_span = 0.02, n_iter = 1)

# rm(list=c('AvRT','Error','J','Length','R','Row','RT','S','ShuffleOrder','LastSubstance','Variation'))

# save(Chromatograms,file="Aligned_Reference_W4_1_Run.RData")




# Post Processing after running the algorithm repeatedly 
# ------------------------------------------------------------------------------------------------------------
# rm(list=ls())
# source('MatrixOperations.R')
# source('MergeRows.R')
# source('ChromaVariation.R')
# load('Aligned_Reference_W4_1_Run.RData')

# Find rows containing no substances --> Mean of rows == NA, remove them

average_rts <- mean_per_row(chromatograms_aligned)

# delete empty rows (if existing)
chromatograms <- lapply(chromatograms_aligned, function(x) {
                                            keep_rows <- which(!is.na(average_rts))
                                            out <- x[keep_rows, ]
                                        })

# still empty rows?
rt_mat1 <- do.call(cbind, lapply(chromatograms, function(x) x$RT))
rowSums(rt_mat1>0)

source("R/merge_rows.R")
average_rts <- mean_per_row(chromatograms)

# min distance here is crucial
chroma_merged <- merge_redundant_rows(chromatograms, average_rts, min_distance=0.02)

average_rts <- mean_per_row(chroma_merged)
rt_mat2 <- do.call(cbind, lapply(chroma_merged, function(x) x$RT))
# save(Chromatograms,file = 'Aligned_Reference_W4_1_Run_Merged_Zeros.RData')

# ----------------------------------------------------------------------------------------------------------------
# # Merge close rows, when only a small proportion of samples contain two peaks,
# # and prevent the merging in the previous run
# 
# Merging <- 'Start'
# while(Merging!='Stop'){
#   # Find and delete Rows that are 
#   # (i) similar in their Retention Time
#   # (ii) redundant i.e. no chromatogram contains substances in both rows
#   # Always start to pick the fist potential candidate of row-pairs to merge
#   # in case merging is not aplicable, take the next position
#   # after merging of one pair, update criteria
#   # if non is ready to merge anymore, stop the merging algorithm
#   
#   # Updating merging criterions
#   average_rts <- mean_per_row(Chromatograms) # Average RTs, after merging rows
#   similar <- SimilarRows(average_rts)    # remaining similarities
#   counter <- 1
#   while (counter!='Stop'){
#     # loop through vector of similar rows, until one shift was done
#     Redundant <- sapply(lapply(Chromatograms, CheckRedundancy,similar[counter]),as.vector) #Check first position
#     Criterion <- IsRedundant(similar = similar[counter],Redundant=Redundant,Criterion = "Proportional")
#     if(Criterion==1){ # only merge if criterion proves redundancy of one of the rows
#       Chromatograms <- lapply(Chromatograms,Merge, ToMerge=similar[counter],Criterion="Area")
#       counter <- 'Stop'
#     } else{
#       counter <- counter+1
#       if (counter>length(similar)){
#         Merging <- 'Stop'
#       }  
#     }
#   }
# }
# 
# save(Chromatograms,file = 'Aligned_Reference_W4_1_Run_Merged_Proportional.RData')

# ---------------------------------------------------------------------------------------------------------
average_rts <- mean_per_row(chroma_merged)

# delete empty rows
del_empty_rows <- function(chromatogram, average_rts){
    chromatogram <- chromatogram[!is.na(average_rts), ]
    chromatogram
}

chromatograms <- lapply(chromatograms, del_empty_rows, average_rts)

# delete blanks
blanks <- c("w4", "w3")

# delete one blank
delete_blank <- function(blank, chromatograms) {
    del_substances <- which(chromatograms[[blank]]$RT > 0)
    chroma_out <- lapply(chromatograms, function(x) x[-del_substances, ])    
}

# delete all blanks
for (i in blanks) {
    chromatograms <- delete_blank(i, chromatograms)
}

# delete single substances
# create matrix with all retention times
rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x$RT))
# find single retention times in rows
single_subs_ind <- which(rowSums(rt_mat > 0) == 1)
# delete substances occuring in just one individual
chromatograms <- lapply(chromatograms, function(x) x[-single_subs_ind, ]) 



# Define the Retention times of each chromatogram by the mean retention time among all 
# rm(list=ls())
# load('Aligned_Reference_W4_1_Run_Merged_Proportional.RData')
# source('ChromaVariation.R')
# source('CleanChromas.R')
# source('MatrixOperations.R')
# source('RelativeAbundance.R')
# average_rts <- mean_per_row(chroma_merged)
# chromatograms <- lapply(chroma_merged,EqualRTs,average_rts) # Change RTs to mean of given row
# 
# Chromatograms <- lapply(Chromatograms,DeleteEmptyRows,average_rts) # Remove empty rows
# Blanks <- BlankPeaks(Chromatograms)
# Chromatograms <- lapply(Chromatograms,DeleteSubstances,Blanks)
# 
# SingleChroma <- SingleSubstance(Chromatograms) # Substances present in just one chromatogram
# Chromatograms <- lapply(Chromatograms,DeleteSubstances,SingleChroma)
# 
# j <- lapply(Chromatograms,ExtractRT)
# 
# 
# save(Chromatograms,file = 'Aligned_Cleaned.RData')
# #load('Aligned_Cleaned.RData')
# 
# Chromatograms <- lapply(Chromatograms,RelativeAbundance)
# save(Chromatograms,file='Preen_Scent.RData')
