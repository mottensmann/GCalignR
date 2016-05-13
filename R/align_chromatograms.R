#'Aligning chromatograms based on retention times
#'
#'@param datafile path to the datafile. It has to be a .txt file. The first row needs to contain
#'  sample names, the second row column names of the corresponding chromatograms. Starting with the
#'  third row chromatograms are included, whereby single samples are concatenated horizontally Each
#'  chromatogram needs to consist of the same number of columns, at least two are required (the
#'  retention time and the area).
#'
#'@param sep  the field separator character. Values on each line of the file are separated by this
#'  character. The default is tab seperated (sep = "\t"). See sep argument in read.table for details.
#'
#'@param rt_name Name of the column holding retention times (i.e. "RT")
#'
#'@param write_output If specified the aligned output is written to a text file. (i.e. write_output = c("RT", "area")
#'       to write the aligned matrices of the retention times and areas to txt files (the names have
#'       to correspond to the column names given in the second row of the datafile).
#'
#'@param rt_cutoff_low lower threshold under which retention times are cutted (i.e. 5)
#'
#'@param rt_cutoff_high Upper threshold above which retention times are cutted (i.e. 35)
#'
#'@param reference A sample(name) to which all other samples are aligned to by means of a
#'  linear shift (i.e. "individual3") The name has to correspond to an individual name given
#'  in the first line of the datafile.
#'
#'@param step1_maxshift Defines a window to search for an optimal linear shift of samples
#'  with respect to the reference. Shifts are evaluated within - step1_maxshift to + step1_maxshift.
#'  Default is 0.05.
#'
#'@param step2_maxshift Defines the allowed deviation of retention times around the mean of
#'  the corresponding row. Default is 0.02.
#'
#'
#' @param blanks character vector of blanks. If specified, all substance found in any of the blanks
#' will be removed from all samples
#'
#'@param step3_maxdiff Defines the minimum difference in retention times among distinct
#'  substances. Substances that do not differ enough, are merged if applicable. Defailt is 0.05.
#'
#'
#'@param blanks Character vector of the names of blanks. If specified, all substance found in any of the blanks
#'  will be removed from all samples (i.e. c("blank1", "blank2")). The names have to correspond
#'  to a name given in the first line of the datafile.
#'
#'@param delete_single_sub logical, determines whether substances that occur in just one sample are
#'  removed or not.
#'
#'
#'@return
#'
#'
#'@references
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@import stringr
#'@import readr
#'
#'@export
#'
#'

align_chromatograms <- function(datafile, sep = "\t", rt_name = NULL, write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
                                step1_maxshift = 0.05, step2_maxshift = 0.02, step3_maxdiff = 0.05, blanks = NULL,
                                del_single_sub = FALSE) {

    if (is.null(rt_name)) stop("specify name of retention time column")
    if (is.null(reference)) stop("specify a reference chromatogram to align the others to")

    # extract names
    ind_names <- readr::read_lines(datafile, n_max = 1) %>%
        stringr::str_split(pattern = sep) %>%
        unlist() %>%
        .[. != ""]

    # extract column names
    col_names <- readr::read_lines(datafile, n_max = 1, skip = 1) %>%
        stringr::str_split(pattern = sep) %>%
        unlist() %>%
        .[. != ""]

    # remove leading and tailing whitespaces
    col_names <- stringr::str_trim(col_names)
    ind_names <- stringr::str_trim(ind_names)

    # extract data
    chroma <- read.table(datafile, skip = 2, sep = sep, stringsAsFactors = F)

    # remove pure NA rows
    chroma <- chroma[!(rowSums(is.na(chroma)) == nrow(chroma)), ]

    # for cara replace commas
    # chroma <- apply(chroma, 2, function(x) x <- str_replace_all(x, ",", "."))

    # transform to numeric
    chroma <-  as.data.frame(apply(chroma, 2, as.numeric))

    # check 1
    if (!((ncol(chroma) / length(col_names)) %% 1) == 0) stop("Number of data columns is not a multiple of the column names provided")
    # check 2
    if (!((ncol(chroma) / length(col_names))  == length(ind_names))) stop("Number of sample names provided does not fit to the number of columns in the data")

    # matrix to list
    chromatograms <- conv_gc_mat_to_list(chroma, ind_names, var_names = col_names)

    #####################
    # save for later use
    ####################
    chroma_raw <- lapply(chromatograms,matrix_append,chromatograms)

    ########################
    # Start of processing
    ########################

    # 1.) cut retention times
    chromatograms <- lapply(chromatograms, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_name)

    # 2.) Linear Transformation of Retentiontimes

    ## thinking about reference: default is chromatogram with most peaks - optional: manual  !Currently it is manual
    chroma_aligned <- linear_transformation(chromatograms, shift=step1_maxshift, step_size=0.01,
                                            error=0, reference = reference, rt_col_name = rt_name)

    # Make List equal in length
    # source("R/matrix_append.R")
    chromatograms <- lapply(chroma_aligned, matrix_append, chroma_aligned)

    ######################
    # save for later use
    ####################
    chroma_linear <- chromatograms

    # source("R/evaluate_chroma.R")
    # Length <- (max(unlist(lapply(chromatograms, function(x) out <- nrow(x))))) # To obtain Rows after run of the algorithm
    # Variation <- mean(var_per_row(chromatograms),na.rm = T)

    # align peaks
    chromatograms_aligned <- align_individual_peaks(chromatograms, error_span = step2_maxshift, n_iter = 1, rt_col_name = rt_name)


    # see whether zero rows are present
    average_rts <- mean_per_row(chromatograms_aligned, rt_col_name = rt_name)

    # delete empty rows (if existing)
    chromatograms <- lapply(chromatograms_aligned, function(x) {
        keep_rows <- which(!is.na(average_rts))
        out <- x[keep_rows, ]
    })

    # calculate average rts again for merging
    average_rts <- mean_per_row(chromatograms, rt_col_name = rt_name)

    # merging step
    # min distance here is crucial --> depends on sample size
    chroma_merged <- merge_redundant_rows(chromatograms, average_rts, min_distance=step3_maxdiff, rt_col_name = rt_name)

    ### just evaluation
    # average_rts <- mean_per_row(chroma_merged)
    # rt_mat2 <- do.call(cbind, lapply(chroma_merged, function(x) x$RT))

    average_rts <- mean_per_row(chroma_merged, rt_col_name = rt_name)

    # delete empty rows again
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
    rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x[[rt_name]]))

    if (del_single_sub) {
        # find single retention times in rows
        single_subs_ind <- which(rowSums(rt_mat > 0) == 1)
        # delete substances occuring in just one individual
        chromatograms <- lapply(chromatograms, function(x) x[-single_subs_ind, ])
    }

    # calculate final retention times
    rt_mat <- do.call(cbind, lapply(chromatograms, function(x) x[[rt_name]]))

#     #######################
#     # Outputs for Heatmaps
#     #######################
#
#     rt_raw <- rt_extract_heatmap(chromatograms = chroma_raw,blanks = blanks,rt_name = rt_name,del_single_sub = del_single_sub)
#     rt_linear <- rt_extract_heatmap(chromatograms = chroma_linear,blanks = blanks,rt_name = rt_name,del_single_sub = del_single_sub)
    rt_raw <- rt_extract_heatmap(chromatograms = chroma_aligned,blanks = blanks,rt_name =rt_name,del_single_sub=del_single_sub)
    rt_linear <- rt_extract_heatmap(chromatograms = chroma_linear,blanks = blanks,rt_name =rt_name,del_single_sub=del_single_sub)
    rt_aligned <- rt_extract_heatmap(chromatograms = chromatograms,blanks = blanks,rt_name =rt_name,del_single_sub=del_single_sub)


    # ============================
    # mean per row without 0
    row_mean <- function(x) {
        if (all(x==0)) 0 else mean(x[x!=0])
    }
    mean_per_row <- apply(rt_mat,1, row_mean)


    # create output matrices for all variables
    output <- lapply(col_names, function(y) as.data.frame(do.call(cbind, lapply(chromatograms, function(x) x[y]))))
    output <- lapply(output, function(x){
                        names(x) <- ind_names
                        x
    })

    output <- lapply(output, function(x){
        x <- cbind(mean_per_row, x)
        x
    })

    output <- lapply(output, function(x){
                        names(x)[1] <- "mean_RT"
                        x
    })
    names(output) <- col_names


    if (!is.null(write_output)){
        write_files <- function(x) {
            write.table(output[[x]], file = paste0(x, ".txt"), sep = "\\t", row.names = FALSE)
        }
        lapply(write_output, write_files)
    }

    #output
    output_algorithm <- list(call=match.call(),
                chroma_aligned = output, rt_raw = rt_raw, rt_linear = rt_linear,
                rt_aligned = rt_aligned)

    class(output_algorithm) <- "GCalign" # name of list
    return(output_algorithm)
}


