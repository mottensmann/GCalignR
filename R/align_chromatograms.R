#' Aligning chromatograms based on retention times
#'
#' @param datafile
#' @param rt_name
#'
#' @return
#'
#'
#' @references
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com) &
#'         Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#' @export
#'


align_chromatograms <- function(datafile, rt_name = "RT", rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
                                step1_maxshift = 0.05, step2_maxshift = 0.02, step3_maxdiff = 0.05, blanks = NULL,
                                del_single_sub = FALSE) {

    # extract names
    ind_names <- readr::read_lines(datafile, n_max = 1) %>%
        stringr::str_split(pattern = "\t") %>%
        unlist() %>%
        .[. != ""]

    # extract column names
    col_names <- readr::read_lines(datafile, n_max = 1, skip = 1) %>%
        stringr::str_split(pattern = "\t") %>%
        unlist() %>%
        .[. != ""]

    # extract data
    chroma <- read.table(datafile, skip = 2, sep = "\t", stringsAsFactors = F)
    # how to deal with commas?
    chroma <- as.data.frame(apply(chroma, 2, function(x) str_replace_all(x, ",", ".")))
    chroma <-  as.data.frame(apply(chroma, 2, as.numeric))

    # check 1
    if (!((ncol(chroma) / length(col_names)) %% 1) == 0) stop("number of data columns is not a multiple of the column names provided")
    # check 2
    if (!((ncol(chroma) / length(col_names))  == length(ind_names))) stop("the number of individual ids provided does not fit to the number of columns in the data")

    # matrix to list
    # source("R/conv_gc_mat_to_list.R")
    # source("R/rename_cols.R")
    chromatograms <- conv_gc_mat_to_list(chroma, ind_names, var_names = col_names)


    # Start of processing --------------------------------------------------------------------------

    # source("R/rt_cutoff.R")
    # 1.) cut retention times below 8
    chromatograms <- lapply(chromatograms, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_name)

    # source("R/linear_transformation.R")
    # 2.) Linear Transformation of Retentiontimes

    ## thinking about reference: default is chromatogram with most peaks - optional: manual
    chroma_aligned <- linear_transformation(chromatograms, shift=step1_maxshift, step_size=0.01, error=0, reference = reference, rt_col_name = rt_name)

    # Make List equal in length
    # source("R/matrix_append.R")
    chromatograms <- lapply(chroma_aligned, matrix_append, chroma_aligned)


    # source("R/evaluate_chroma.R")
    Length <- (max(unlist(lapply(chromatograms, function(x) out <- nrow(x))))) # To obtain Rows after run of the algorithm
    Variation <- mean(var_per_row(chromatograms),na.rm = T)

    #rm(list=c("Chroma_aligned","chroma","AdjustRetentionTime","AlignPeaks","BestShift",
    #          "RetentionCutoff","SharedPeaks","PeakShift"))

    # source("R/align_individual_peaks.R")
    # source("R/shift_rows.R")
    chromatograms_aligned <- align_individual_peaks(chromatograms, error_span = step2_maxshift, n_iter = 1)


    average_rts <- mean_per_row(chromatograms_aligned)

    # delete empty rows (if existing)
    chromatograms <- lapply(chromatograms_aligned, function(x) {
        keep_rows <- which(!is.na(average_rts))
        out <- x[keep_rows, ]
    })

    # still empty rows?
    rt_mat1 <- do.call(cbind, lapply(chromatograms, function(x) x$RT))
    rowSums(rt_mat1>0)

    # source("R/merge_rows.R")
    average_rts <- mean_per_row(chromatograms)

    # min distance here is crucial --> depends on sample size
    chroma_merged <- merge_redundant_rows(chromatograms, average_rts, min_distance=step3_maxdiff)

    average_rts <- mean_per_row(chroma_merged)
    rt_mat2 <- do.call(cbind, lapply(chroma_merged, function(x) x$RT))


    average_rts <- mean_per_row(chroma_merged)

    # delete empty rows
    del_empty_rows <- function(chromatogram, average_rts){
        chromatogram <- chromatogram[!is.na(average_rts), ]
        chromatogram
    }
    chromatograms <- lapply(chromatograms, del_empty_rows, average_rts)


    # delete blanks
    if (!is.null(blanks)) {
        # delete one blank
        delete_blank <- function(blank, chromatograms) {
            del_substances <- which(chromatograms[[blank]]$RT > 0)
            chroma_out <- lapply(chromatograms, function(x) x[-del_substances, ])
        }
        # delete all blanks
        for (i in blanks) {
            chromatograms <- delete_blank(i, chromatograms)
        }
    }


    # delete single substances
    # create matrix with all retention times
    rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x$RT))

    if (del_single_sub) {
        # find single retention times in rows
        single_subs_ind <- which(rowSums(rt_mat > 0) == 1)
        # delete substances occuring in just one individual
        chromatograms <- lapply(chromatograms, function(x) x[-single_subs_ind, ])
    }

    # calculate final retention times
    rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x$RT))

    # mean per row without 0
    row_mean <- function(x) {
        if (all(x==0)) 0 else mean(x[x!=0])
    }
    mean_per_row <- apply(rt_mat,1, row_mean)

    # create abundance matrix
    area_df <- as.data.frame(do.call(cbind, lapply(chromatograms, function(x) x$Area)))

    # rownames
    row.names(area_df) <- mean_per_row

}


