#' Merge two rows into one
#'
#' @description
#' Merges the content of two rows in such a way that content is retained. This means that either the row
#' is selected that contains in case that the second is empty (i.e just zeros). If both rows have a content
#' the row with the larger concentration (i.e peak area or peak height) is retained in \code{gc_peak_df}
#'
#' @inheritParams align_chromatograms
#'
#' @inheritParams matrix_append
#'
#' @param to_merge
#' integer vector indicating the index of the last row of a pair to merge.
#'
#' @param criterion
#' character string either "zero" or "area". "zero" is used to delete the row containing zeros,
#' "area" provokes the selection of the larger of two peaks, defined by the concentration of the peaks.
#' "area" is only useful when \code{criterion} is "proportional".
#'
#'
#' @return
#' aligned chromatograms
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
#'
merge_rows <- function(gc_peak_df, to_merge, criterion="zero", rt_col_name){
    # Check always the row containing just zeros, in case of zeros in both, just delete one of them
    # To Merge == Last row of a similar pair
    Row1 <- to_merge-1 # Previous Row
    Row2 <- to_merge # Current Row
    R1 <- gc_peak_df[Row1, rt_col_name]
    R2 <- gc_peak_df[Row2, rt_col_name]
    if (criterion=="zero"){
       if (Row1 > 1 & Row2 < nrow(gc_peak_df)){ # Avoid taking first and last rows
            if (R1 == 0){
                #  Delete Row1
                gc_peak_df <- rbind(gc_peak_df[1:(Row1-1), ], gc_peak_df[Row2:nrow(gc_peak_df), ])
            } else if (R2 == 0){
                # Delete Row2
                gc_peak_df <- rbind(gc_peak_df[1:Row1,],gc_peak_df[(Row2+1):nrow(gc_peak_df), ])
            }
       }

      if (Row2 == nrow(gc_peak_df)){
          if (R1 == 0) {
            gc_peak_df <- rbind(gc_peak_df[1:(Row1-1), ], gc_peak_df[Row2:nrow(gc_peak_df), ])
        } else if (R2 == 0){
            gc_peak_df <- gc_peak_df[1:Row1, ]
        }
      }
    }

    if(criterion=="area"){ # Take the larger of two peaks, defined by the are of the peaks
        if (gc_peak_df$Area[Row1] >= gc_peak_df$Area[Row2]){
            gc_peak_df <- rbind(gc_peak_df[1:Row1,] ,gc_peak_dfs[(Row2+1):nrow(gc_peak_df), ])
        } else if (gc_peak_df$Area[Row1] < gc_peak_df$Area[Row2]) {
            gc_peak_df <- rbind(gc_peak_df[1:(Row1-1),],gc_peak_df[Row2:nrow(gc_peak_df), ])

        }
    }
    gc_peak_df
}
