#' calculates mean of retention times up to a certain sample (not including 0s)
#' 
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param samples indices of samples up to sample of interest (1:sample-1)
#' @param retention_row current retention time row to be compared
#' 
#' @return 
#' mean of rts 
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      

mean_of_samples <- function(chromatograms, samples, retention_row){
    rts <- unlist(lapply(chromatograms[samples], function(x) x$RT[retention_row]))
    mean_rt <- mean(rts[!(rts == 0)], na.rm = TRUE)
    ## round ?
    
}

#' calculates mean of retention times up to a certain sample (not including 0s)
#' 
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param samples indices of samples up to sample of interest (1:sample-1)
#' @param retention_row current retention time row to be compared
#' 
#' @return 
#' var of rts 
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      
var_of_samples <- function(chromatograms, samples, retention_row){
    # Estimate the Variation in Retention Times within Rows  
    # NA indicates that only one Substance exists  
    rts <- unlist(lapply(chromatograms[samples], function(x) x$RT[retention_row]))
    mean_rt <- var(rts[!(rts == 0)], na.rm = TRUE)
}


#' calculates mean retention time per row
#' 
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @return 
#' vector with mean retention times per row in chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      
# Evaluate Algorithm performance
mean_per_row = function(chromatograms){
    n_substance <- nrow(chromatograms[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance, 
                    function(x) mean_of_samples(chromatograms, 1:length(chromatograms), x)))
    out

}


#' calculates mean retention time per row
#' 
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @return 
#' vector with variance of retention times per row in chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'      
var_per_row = function(chromatograms){
    n_substance <- nrow(chromatograms[[1]]) # all_chromatograms have equal number of rows
    out <- unlist(lapply(1:n_substance, 
        function(x) var_of_samples(chromatograms, 1:length(chromatograms), x)))
    out
    
}







# 
# SimilarRows = function(AverageRTs,MinDistance=0.05){
#     # Estimate Degree of Similarity between subsequent rows by comparison of mean
#     # retention times of rows.
#     # Similarity is evaluated at the level of MinDistance, given in seconds
#     # The output gives the the position of a row, that is similar to the previous
#     # i.e. row 5 is similar to 4
#     
#     Difference <- rep(NA, (length(AverageRTs)-1)) # Estimate Difference between adjacent rows
#     for (i in 2:length(AverageRTs)){
#         Difference[i] <- AverageRTs[i]-AverageRTs[i-1]  
#     }
#     similar <- which(Difference<=MinDistance) # Which rows differ less the MinDistance ?
#     similar
# }
# 
# CheckRedundancy = function(Chromatograms,similar){
#     # If only one of two neighbouring rows contain a substance
#     # they are redundant, coded by a One
#     Row1 <- Chromatograms$RT[similar-1] # Extract previous row
#     Row2 <- Chromatograms$RT[similar] # Extract current row
#     Redundant <- 0
#     if (Row1==0 | Row2==0){
#         Redundant <- 1
#     }
#     Redundant
# }
# 
# IsRedundant = function(similar, Redundant, Criterion="Strict"){
#     # Indicates by a binary output variable (1/0) if rows should be merged
#     # Methods: Strict: A single sample with two peaks prevents merging
#     #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
#     ToMerge <- 0
#     if (Criterion=="Strict"){
#         if(sum(Redundant)/length(Redundant)==1){
#             ToMerge <- 1
#         }
#     } else if (Criterion=="Proportional"){
#         if(sum(Redundant)/length(Redundant)>=0.95){
#             ToMerge <- 1
#         }
#     }
#     ToMerge
# }
# 
# MeanOfSamples = function(Data,Objects, Row){
#     # Calculate mean retention times within a row for a number of objects
#     # Data is list of Chromatograms
#     # Objects is a vector specifying the objects on which the calculation of the mean is based
#     # Row indicates the row to look at
#     
#     AvRT <- numeric(0)
#     for (N in 1:length(Objects)){
#         AvRT <- c(AvRT,Data[[Objects[N]]]$RT[Row])
#     }  
#     AvRT <- mean(AvRT[AvRT>0])       # Neglect Zeros (i.e. do not incorporate samples without the substance in the Row)
#     AvRT <- round(AvRT,digits = 2)
# }