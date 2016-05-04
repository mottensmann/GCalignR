#' merges redundant rows
#'
#' @param chromatograms \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param average_rts average retention times across samples per row
#' @param min_distance difference between the mean retention time of two rows of the chromatograms
#'        to be considered for merging if no individual has substances in both rows.
#'
#' @return
#' chromatograms with merged rows
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
merge_redundant_rows <- function(chromatograms, average_rts,  min_distance=0.05){

    Merging <- 'Start'
    while(Merging != 'Stop'){

        # Find and delete Rows that are
        # (i) similar in their Retention Time
        # (ii) redundant i.e. no chromatogram contains substances in both rows
        # Always start to pick the fist potential candidate of row-pairs to merge
        # in case merging is not aplicable, take the next position
        # after merging of one pair, update criteria
        # if non is ready to merge anymore, stop the merging algorithm

        # Updating merging criterions
        average_rts <- mean_per_row(chromatograms) # Average RTs, after merging rows
        similar <- similar_rows(average_rts, min_distance)    # remaining similarities
        # print(similar)

       #  if (sum(similar > nrow(chromatograms[[1]])) > 0) similar <- similar[!(similar > nrow(chromatograms[[1]]))]
        counter <- 1
        while (counter!='Stop'){
            # loop through vector of similar rows, until one shift was done
            # print(paste0("counter=", counter))

            # check if any chromatogram contains substances in both rows
            if (length(similar) == 0){
                Merging <- "Stop"
                break
            }
            redundant <- sapply(lapply(chromatograms, check_redundancy, similar[counter]),as.vector) #Check first position

            # checks whether all individuals just have substances in one of the rows ("Strict")
            # checks whether at least 95% of individuals just have substances in one of the rows ("Proportional")
            criterion <- is_redundant(redundant = redundant, criterion="strict")

            if (criterion == 1){ # only merge if criterion proves redundancy of one of the rows
                chromatograms <- lapply(chromatograms, merge_rows, to_merge = similar[counter], criterion="zero")

                # check2 <- similar[counter]
                counter <- 'Stop'

               # check <- do.call(cbind, lapply(chromatograms, function(x) x$RT))

            } else if  (criterion == 0) {
                counter <- counter+1

                if (counter>length(similar)){
                    Merging <- 'Stop'
                    counter <- 'Stop'
                }
            }
        }
    }

    chromatograms
}


#' Estimate degree of similarity between subsequent rows by comparison of mean retention times of rows.
#'
#' @param average_rts average retention times per row across individuals
#' @param  min_distance default 0.05,
#' @details Similarity is evaluated at the level of min_distance, given in seconds
#' @return
#' position of a row, that is similar to the previous, i.e. row 5 is similar to 4
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
#'
similar_rows <- function(average_rts, min_distance=0.05){

    difference <- rep(NA, (length(average_rts)-1)) # Estimate Difference between adjacent rows

    for (i in 2:length(average_rts)){
        difference[i] <- average_rts[i]-average_rts[i-1]
    }
    similar <- which(difference <= min_distance) # Which rows differ less the MinDistance ?
    similar
}

#' If only one of two neighbouring rows contain a substancethey are redundant, coded by a One
#'
#' @param chromatograms
#' @param  similar
#'
#' @return
#' if similar rows are redundant, 1, else 0
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'

check_redundancy <- function(chromatogram, similar){
    # If only one of two neighbouring rows contain a substance
    # they are redundant, coded by a One
    Row1 <- chromatogram$RT[similar-1] # Extract previous row
    Row2 <- chromatogram$RT[similar] # Extract current row
    Redundant <- 0
    if (Row1==0 | Row2==0){
        Redundant <- 1
    }
    Redundant
}

#' Indicates by a binary output variable (1/0) if rows should be merged
#'
#' @param chromatograms
#' @param  similar
#'
#' @details  Methods: "strict": A single sample with two peaks prevents merging
#           "proportional": Merging is acceptabel if only 5 % of samples show two peaks
#' @return
#' if similar rows are redundant, 1, else 0
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
#'
is_redundant <- function(similar, redundant, criterion="strict"){
    # Indicates by a binary output variable (1/0) if rows should be merged
    # Methods: Strict: A single sample with two peaks prevents merging
    #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
    ToMerge <- 0
    if (criterion == "strict"){
        if(sum(redundant)/length(redundant) == 1){
            ToMerge <- 1
        }
    } else if (criterion == "proportional"){
        if(sum(redundant)/length(redundant)>=0.95){
            ToMerge <- 1
        }
    }
    ToMerge
}

#' algorithm to adjust unreliability of chromatograms on an individual peak basis
#'
#' @param chromatogram \code{data.frame} containing GC data (retention times, peak area, peak hight etc) for
#'   one individual in adjacent columns. The first column for all individuals has to be the retention
#'   time, retention time has to be named RT.
#' @param to_merge Last row of a similar pair
#' @param criterion "Zero" - delete a row containing zeros, or
#'        "area": take the larger of two peaks, defined by the are of the peaks. "area" is just
#'        useful when redundancy criterion was "proportional" instead of the default "strict".
#' @return
#' aligned chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
merge_rows <- function(chromatogram, to_merge, criterion="zero"){
    # Check always the row containing just zeros, in case of zeros in both, just delete one of them
    # To Merge == Last row of a similar pair
    Row1 <- to_merge-1 # Previous Row
    Row2 <- to_merge # Current Row
    R1 <- chromatogram$RT[Row1]
    R2 <- chromatogram$RT[Row2]
    if (criterion=="zero"){
       if (Row1 > 1 & Row2 < nrow(chromatogram)){ # Avoid taking first and last rows
            if (R1 == 0){
                #  Delete Row1
                chromatogram <- rbind(chromatogram[1:(Row1-1), ], chromatogram[Row2:nrow(chromatogram), ])
            } else if (R2 == 0){
                # Delete Row2
                chromatogram <- rbind(chromatogram[1:Row1,],chromatogram[(Row2+1):nrow(chromatogram), ])
            }
       }

      if (Row2 == nrow(chromatogram)){
          if (R1 == 0) {
            chromatogram <- rbind(chromatogram[1:(Row1-1), ], chromatogram[Row2:nrow(chromatogram), ])
        } else if (R2 == 0){
            chromatogram <- chromatogram[1:Row1, ]
        }
      }
    }

    if(criterion=="area"){ # Take the larger of two peaks, defined by the are of the peaks
        if (chromatogram$Area[Row1] >= chromatogram$Area[Row2]){
            chromatogram <- rbind(chromatogram[1:Row1,] ,chromatograms[(Row2+1):nrow(chromatogram), ])
        } else if (chromatogram$Area[Row1] < chromatogram$Area[Row2]) {
            chromatogram <- rbind(chromatogram[1:(Row1-1),],chromatogram[Row2:nrow(chromatogram), ])

        }
    }
    chromatogram
}
