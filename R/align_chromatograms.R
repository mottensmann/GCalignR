#' Aligning peaks based on retention times
#'
#'@description
#' This is the core function of \code{\link{GCalignR}} to align peak data. The input data is a peak list. Read through the documentation below and take a look at the vignettes for a thorough introduction. Three parameters \code{max_linear_shift}, \code{max_diff_peak2mean} and \code{min_diff_peak2peak} are required as well as the column name of the peak retention time variable \code{rt_col_name}. Arguments are described among optional parameters below.
#'
#'@details
#' This function aligns and matches homologous peaks across samples using a three-step algorithm based on user-defined parameters that are explained in the next section. In brief: \strong{(1)} A full alignment of peak retention times is conducted to correct for systematic linear drift of retention times among homologous peaks from run to run. Thereby a coarse alignment is achieved that maximises the similarity of retention times across homologous peaks prior to a \strong{(2)} partial alignment and matching of peaks. This and the next step in the alignment is based on a retention time matrix that contains all samples as columns and peak retention times in rows. The goal is to match homologous peaks within the same row that represents a chemical substance. Here, peaks are recognised as homologous when the retention time matches within a user-defined error span (see \code{max_diff_peak2mean}) that approximates the expected retention time uncertainty. Here, the position of every peak in the matrix is evaluated in comparison to similar peaks across the complete dataset and homologous peaks are gradually grouped together row by row. After all peaks were matched, a \strong{(3)} adjacent rows of similar retention time are scanned to detect redundancies. A pair of rows is identified as redundant and merged if mean retention times are closer than specified by \code{min_diff_peak2peak} and single samples only contain peak in one but not both rows. Therefore, this step allows to match peaks that are associated with higher variance than allowed during the previous step. Several optional processing steps are available, ranging from the removal of peaks representing contaminations (requires to include blanks as a control) to the removal of uninformative peaks that are present in just one sample (so called singletons).
#'
#'@inheritParams check_input
#'
#'
#'@param data
#' Dataset containing peaks that need to be aligned and matched. For every peak a arbitrary number of numerical variables can be included (e.g. peak height, peak area) in addition to the mandatory retention time. The standard format is a tab-delimited text file according to the following layout: (1) The first row contains sample names, the (2) second row column names of the corresponding peak lists. Starting with the third row, peak lists are included for every sample that needs to be incorporated in the dataset. Here, a peak list contains data for individual peaks in rows, whereas columns specify variables in the order given in the second row of the text file. Peak lists of individual samples are concatenated horizontally and need to be of the same width (i.e. the same number of columns in consistent order). Alternatively, the input may be a list of data frames. Each data frame contains the peak data for a single individual. Variables (i.e.columns) are named consistently across data frames. The names of elements in the list are used as sample identifiers. Cells may be filled with numeric or integer values but no factors or characters are allowed. NA and 0 may be used to indicate empty rows.
#'
#'@param sep
#' The field separator character. The default is tab separated (\code{sep = '\\t'}).
#' See the "sep" argument in \code{\link[utils]{read.table}} for details.
#'
#'@param rt_col_name
#' A character giving the name of the column containing the retention times. The decimal separator needs to be a point.
#'
#'@param write_output
#' A character vector of variable names. For each variable a text file containing the aligned dataset is written to the working directory. Vector elements need to correspond to column names of data.
#'
#'@param rt_cutoff_low
#' A numeric value giving a retention time threshold. Peaks with retention time below the threshold are removed in a preprocessing step.
#'
#'@param rt_cutoff_high
#' A numeric value giving a retention time threshold. Peaks with retention time above the threshold are removed in a preprocessing step.
#'
#'@param reference
#' A character giving the name of sample from the dataset. By default, a sample is automatically selected from the dataset using the function \code{\link{choose_optimal_reference}}. The reference is used for the full alignment of peak lists by linear transformation.
#'
#'@param max_linear_shift
#' Numeric value giving the window size considered in the full alignment. Usually, the amplitude of linear drift is small in typical GC-FID datasets. Therefore, the default value of 0.05 minutes is adequate for most datasets. Increase this value if the drift amplitude is larger.
#'
#'@param max_diff_peak2mean
#' Numeric value defining the allowed deviation of the retention time of a given peak from the mean of the corresponding row (i.e. scored substance). This parameter reflects the retention time range in which peaks across samples are still matched as homologous peaks (i.e. substance). Peaks with retention times exceeding the threshold are sorted into a different row.
#'
#'@param min_diff_peak2peak
#' Numeric value defining the expected minimum difference in retention times among homologous peaks (i.e. substance). Rows that differ less in the mean retention time, are therefore merged if every sample contains either one or none of the respective compounds. This parameter is a major determinant in the classification of distinct peaks. Therefore careful consideration is required to adjust this setting to your needs (e.g. the resolution of your gas-chromatography pipeline). Large values may cause to merge truly different substances with similar retention times, if those are not simultaneously occurring within at least one individual, which might occur by chance for small sample sizes. By default set to 0.2 minutes.
#'
#'@param blanks
#' Character vector of names of negative controls. Substances found in any of the blanks will be removed from the aligned dataset, before the blanks are deleted from the aligned data as well. This is an optional filtering step.
#'
#'@param delete_single_peak
#' Boolean, determining whether substances that occur in just one sample are removed or not.
#'
#' @param remove_empty
#' Boolean, allows to remove samples which lack any peak after the alignment finished. By default FALSE
#'
#' @param permute
#' Boolean, by default a random permutation of samples is conducted prior for each row-wise alignment step. Setting this parameter to FALSE causes alignment of the dataset as it is.
#'
#' order of samples is constantly randomised during the alignment. Allows to prevent this behaviour for maximal repeatability if needed.
#'
#'@return
#' Returns an object of class "GCalign" that is a a list containing several objects that are listed below. Note, that the objects "heatmap_input" and "Logfile" are best inspected by calling the provided functions \code{gc_heatmap} and \code{print}.
#'
#'
#' \item{aligned}{Aligned Gas Chromatography peak data subdivided into individual data frames for every variable. Samples are represented by columns, rows specify homologous peaks. The first column of every data frame is comprised of the mean retention time of the respective peak (i.e. row). Retention times of samples resemble the values of the raw data. Internally, linear adjustments are considered where appropriate}
#' \item{heatmap_input}{Used internally to create heatmaps of the aligned data}
#' \item{Logfile}{A protocol of the alignment process.}
#' \item{input_list}{Input data in form of a list of data frames.}
#' \item{aligned_list}{Aligned data in form of a list of data frames.}
#' \item{input_matrix}{List of matrices. Each matrix contains the input data for a variable}
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'
#'@examples
#' ## Load example dataset
#' data("peak_data")
#' ## Subset for faster processing
#' peak_data <- peak_data[1:3]
#' peak_data <- lapply(peak_data, function(x) x[1:50,])
#' ## align data with default settings
#' out <- align_chromatograms(peak_data, rt_col_name = "time")
#'
#'@export
#'
align_chromatograms <- function(data, sep = "\t", rt_col_name = NULL,
                                write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
                                max_linear_shift = 0.02, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.08, blanks = NULL, delete_single_peak = FALSE, remove_empty = FALSE, permute = TRUE, ...) {

    opt <- list(...)

    # Print start
    cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%Y-%m-%d %H:%M:%S")),'\n\n'))

    # Iteration have been deleted as function paramter.
    iterations = 1

    ### 1. Checks preparations
    ### ======================

    # 1.1 Stop execution if mandatory checks are not passed
    if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
    x <- check_input(data,sep,write_output = write_output,blank = blanks,reference = reference,rt_col_name = rt_col_name, message = FALSE)
    if (x != TRUE) stop("Processing not possible: check warnings below and change accordingly in order to proceed")

    # 1.2 Create a "Logbook" to record alignment steps and parameters
    Logbook <- list()
    Logbook[["Date"]]["Start"] <- as.character(strftime(Sys.time()))

    ### 1.3. Check that parameter combinations are sensible
    ### ===================================================
    # if (min_diff_peak2peak <= max_diff_peak2mean) {
    #     cat("It is not advisable to set 'min_diff_peak2peak' to the same or a smaller value than 'max_diff_peak2mean'. See the vignettes for further information.\n")
    #     askYesNo <- function() {
    #         x <- readline("Do you want to quit the alignment now? [Yes/No] ")
    #         while (!(x %in% c("Yes", "No"))) x <- readline("Type 'Yes' or 'No' ")
    #         if (x == "Yes") stop("Processing interrupted on user request")
    #     }
    #     askYesNo()
    # }

    ### 2. Load Data
    ### ============
    if (is.character(data)) { # txt file
        # Get Sample Names
        ind_names <- readr::read_lines(data, n_max = 1)
        ind_names <- unlist(stringr::str_split(string = ind_names,pattern = sep))
        ind_names <- ind_names[ind_names != ""]
        # Get Variable Names
        col_names <- readr::read_lines(data, n_max = 1, skip = 1)
        col_names <- unlist(stringr::str_split(string = col_names,pattern = sep))
        col_names <- col_names[col_names != ""]
        col_names <- stringr::str_trim(col_names)
        ind_names <- stringr::str_trim(ind_names)
        # Get Data
        gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F, fill = T)
        # Remove just NA-rows
        gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ]
        # Remove empty rows
        gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))]
        # Transform to data frame
        gc_data <-  as.data.frame(gc_data)
        # gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))
        # convert to list
        gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)
        # convert retention times to numeric
        fx <- function(x,rt_col_name) {
            x[[rt_col_name]] <- as.numeric(x[[rt_col_name]])
            return(x)
        }
        gc_peak_list <- lapply(gc_peak_list,FUN = fx,rt_col_name = rt_col_name)
        ## Remove purely-zero rows in samples (added 03.07)
        gc_peak_list <- lapply(gc_peak_list, FUN = function(x) {
            z <- which(x[[rt_col_name]] == 0)
            if (length(z) > 0) x <- x[-z,]
            return(x)
        })


    } else if (is.list(data)) { # data is in a list
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        col_names <- names(data[[1]])
        ind_names <- names(data)
        gc_peak_list <- lapply(data,matrix_append,gc_peak_list = data, val = "NA") # same dimensions of dfs

        # # force to numeric
        # temp <- do.call("cbind", data)
        # if (!(any(apply(temp, 2, class) %in% c("numeric","integer")))) {
        #     na_1 <- length(which(is.na(temp)))
        #     data <- lapply(data, function(x) as.data.frame(apply(x, 2, as.numeric)))
        #     na_2 <- length(which(is.na(do.call("cbind", data))))
        #     if (na_2 > na_1) warning("All columns need to contain only numericals or integers. NAs introduced by coercion")
        #     }

    } # end load data

    ## check that each peak consistently is ordered by increasing
    ## peak retention times
    ordered.input <- sapply(gc_peak_list, function(x) {
        any(diff(order(x[[rt_col_name]])) != 1)
    })
    if (any(ordered.input == TRUE)) {
        gc_peak_list <- lapply(gc_peak_list, function(x) {
            ## sort by increasing retention times
            x[order(x[[rt_col_name]], decreasing = F),]
        })
    }

    ## subset samples
    if ("samples" %in% names(opt)) gc_peak_list[opt[["samples"]]]

    # save gc_peak_list for documentation purposes
    input_list <- gc_peak_list

    # Write some information about the input data to the Logfile
    cat(paste0('Data for ',as.character(length(ind_names)),' samples loaded.'))
    Logbook[["Input"]]["Samples"] <- length(ind_names)

    Logbook[["Input"]]["Range"] <- paste((range(peak_lister(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name))),collapse = "-")
    Logbook[["Input"]]["File"] <- as.character(as.character(match.call()["data"]))
    Logbook[["Input"]]["Retention_Time"] <- rt_col_name
    # Logbook[["Input"]]["Concentration"] <- conc_col_name
    Logbook[["Input"]]["Peaks"] <- peak_counter(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name)
    # Only created if blanks!=NULL
    if (!is.null(blanks)) Logbook[["Input"]][["Blanks"]] <- paste(blanks,collapse = "; ")

    ### 3. Start processing
    ### ===================
    # If the data only contains retention times (i.e. one colum) a dummy variable is added to ensure full compatibility with the data frame based manipulations
    if (ncol(gc_peak_list[[1]]) == 1) {
        gc_peak_list <- lapply(X = gc_peak_list, dummy_col)
        col_names <- c(col_names, "GCalignR_Dummy")
    }
    # 3.1. Cut retention times
    if (!is.null(rt_cutoff_low) | !is.null(rt_cutoff_high)) {
        cat("\nApply retention time cut-offs ... ")
    }
    gc_peak_list <- lapply(gc_peak_list, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_col_name)
    gc_peak_list_raw <- lapply(gc_peak_list, matrix_append, gc_peak_list)
    if (!is.null(rt_cutoff_low) | !is.null(rt_cutoff_high)) {
        cat("done")
    }
    # Write to Logbook, distinguish cases
    if (!is.null(rt_cutoff_low) & is.null(rt_cutoff_high)) {
        Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
        Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
    } else if (!is.null(rt_cutoff_high) & is.null(rt_cutoff_low)) {
        Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
        Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
    } else if (!is.null(rt_cutoff_high) & !is.null(rt_cutoff_low)) {
        Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
        Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
    } else if (is.null(rt_cutoff_high) & is.null(rt_cutoff_low)) {
        Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
        Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
    }

    # 3.2 Linear Transformation of peak retention times

    # Revision 24-04-2017:
    # Rounding is eliminated. Calculations are based on rounded values instead.

    # Round retention times to two decimals
    #round_rt <- function(gc_peak_df,rt_col_name = rt_col_name) {
    #gc_peak_df[rt_col_name] <- round(gc_peak_df[rt_col_name],digits = 2)
    #return(gc_peak_df)
    #}
    # gc_peak_list <- lapply(X = gc_peak_list,FUN = round_rt)

    if (!is.null(reference)) {
        cat(paste0('\nStart correcting linear shifts with ',"\"",as.character(reference),"\"",' as a reference ...\n'))
        gc_peak_list_linear <- linear_transformation(gc_peak_list = gc_peak_list, max_linear_shift = max_linear_shift, step_size = 0.01, reference = reference, rt_col_name = rt_col_name, Logbook = Logbook)
        Logbook <- gc_peak_list_linear[["Logbook"]]
        gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
        gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
        gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames, col_names)
    } else if (is.null(reference)) {
        # If no reference was specified by the user, the reference is determined, such
        # that the sample with the highest avarage similarity to all other samples is used.
        cat("\nNo reference was specified. Hence, a reference will be selected automatically ...\n ")
        best.ref <- choose_optimal_reference(data = gc_peak_list, rt_col_name = rt_col_name)
        # set the reference
        reference <- best.ref[["sample"]]
        cat("\n")
        text <- paste0("'",reference,"' was selected on the basis of highest average similarity to all samples (score = ",round(best.ref[["score"]],2),").\n")
        cat(stringr::str_wrap(paste(text,collapse = ""),width = 100,exdent = 0,indent = 0))
        #    }# new end
        cat(paste0('\nStart correcting linear shifts with ',"\"",as.character(reference),"\"",' as a reference ...\n'))
        gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift = max_linear_shift,
                                                     step_size = 0.01, reference = reference, rt_col_name = rt_col_name, Logbook = Logbook)
        Logbook <- gc_peak_list_linear[["Logbook"]]
        gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
        gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
        gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames,col_names)
    }# old end

    Logbook[["Input"]]["Reference"] <- reference
    # why is there a cat?
    #cat(" done\n")
    # equalise chromatograms sizes
    gc_peak_list_linear <- lapply(gc_peak_list_linear, matrix_append, gc_peak_list_linear)

    # 3.3 Align peaks
    cat("\n")
    cat(c('Start aligning peaks ... ','this might take a while!\n'))

    # create corresponding lists
    gc_peak_list_aligned <- gc_peak_list_linear
    no_peaks <- matrix(NA,nrow = iterations,ncol = 1)
    merged_peaks <- matrix(NA, nrow = iterations,ncol = 1)

    # Iterations are currently limited to one loop
    for (R in 1:iterations) {

        if (max_diff_peak2mean > 0) {
            gc_peak_list_aligned <- align_peaks(gc_peak_list_aligned, max_diff_peak2mean = max_diff_peak2mean,
                                                iterations = iterations, rt_col_name = rt_col_name,R = R,
                                                permute = permute)
        } else if (max_diff_peak2mean == 0) {
            gc_peak_list_aligned <- align_peaks_fast(gc_peak_list_aligned, rt_col_name = rt_col_name)
        }

        # mean rt per row
        average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
        # remove empty rows
        gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) {
            keep_rows <- which(!is.na(average_rts))
            out <- x[keep_rows, ]
        })

        # mean retention time per row
        average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
        # Number of peaks before merging
        no_peaks[R] <- nrow(gc_peak_list_aligned[[1]])

        # 3.4 Merge rows
        cat("\nMerge redundant rows ...\n ")
        gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned,
                                                      min_diff_peak2peak = min_diff_peak2peak, rt_col_name = rt_col_name)
        # estimate Number of merged peaks
        merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])
        # Number of peaks after merging
        no_peaks[R] <- nrow(gc_peak_list_aligned[[1]])

        Logbook[["Aligned"]]["total"] <- no_peaks[R]

        #cat('done')
        average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
        # delete empty rows again
        gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)
    }#end iterative loop of aligning & merging

    #cat(paste0('\n','Alignment completed'),'\n\n')

    ### 4 Cleaning chromatograms
    ### ========================
    gc_peak_list_aligned <- gc_peak_list_aligned[match(names(gc_peak_list_raw),names(gc_peak_list_aligned))]

    # delete peaks present in blanks, then remove the blank itself
    if (!is.null(blanks)) {
        cat("Remove contaminations and remove blanks ... ")
        # delete peaks present in blanks

        # Number of Peaks including Blanks
        N <- nrow(gc_peak_list_aligned[[1]])
        # delete all blanks and peaks found in those

        gc_peak_list_aligned <- delete_blank(blanks = blanks,
                                             gc_peak_list_aligned = gc_peak_list_aligned,
                                             rt_col_name = rt_col_name)

        N <- N -  nrow(gc_peak_list_aligned[[1]])
        Logbook[["Filtering"]]["Blank_Peaks"] <- N
        Logbook[["Aligned"]]["blanks"] <- N
        cat('done\n')
    }

    # Create a retention time matrix
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))
    # To estimate the number of deleted peaks
    singular_peaks <- nrow(rt_mat)
    # find single retention times in rows
    if (delete_single_peak) {
        cat("remove single peaks ... ")
        single_subs_ind <- which(rowSums(rt_mat > 0) == 1)
        # delete substances occurring in just one individual
        if (length(single_subs_ind) > 0) {
            gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) x[-single_subs_ind, ])
        }
        cat(paste(as.character(length(single_subs_ind)),'have been removed\n'))
        Logbook[["Aligned"]]["singular"] <- length(single_subs_ind)
    }

    ## added 2019-11-05
    ## Remove samples which are empty
    if (remove_empty) {
        cat("remove empty samples ...")
        ## check if sum of retention times is zero for entire sample
        check_empty <- lapply(gc_peak_list_aligned, function(x) sum(x[[rt_col_name]]))
        empty_sample_names <- which(check_empty == 0) %>% names
        # remove blanks
        gc_peak_list_aligned[empty_sample_names] <- NULL

        if (length(empty_sample_names) == 0) cat('no sample was removed\n')
        if (length(empty_sample_names) == 1) cat(paste(empty_sample_names), 'was removed\n')
        if (length(empty_sample_names) > 1) cat(paste(empty_sample_names), 'were removed\n')
    }

    # Final RT Matrix
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))
    # how many were deleted, if any
    singular_peaks <- singular_peaks - nrow(rt_mat)
    Logbook[["Aligned"]]["retained"] <- nrow(gc_peak_list_aligned[[1]])

    ### 5 Create Output for Heatmaps
    ### ============================

    # Initial input
    rt_raw <- rt_extract(gc_peak_list = gc_peak_list_raw,rt_col_name = rt_col_name)
    # after linear shifts were corrected
    rt_linear <- rt_extract(gc_peak_list = gc_peak_list_linear,rt_col_name = rt_col_name)
    # final output
    rt_aligned <- rt_extract(gc_peak_list = gc_peak_list_aligned,rt_col_name = rt_col_name)

    ### 6: Create output matrices for Variables
    ### =======================================

    # Aligned data
    row_mean <- function(x) {if (all(x == 0)) 0 else mean(x[x != 0])}
    mean_per_row <- apply(rt_mat,1, row_mean)
    output <- lapply(col_names, function(y) as.data.frame(do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[y]))))
    output <- lapply(output, function(x){
        names(x) <- names(gc_peak_list_aligned)
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

    # Input
    rt_mat <- do.call(cbind, lapply(gc_peak_list_raw, function(x) x[[rt_col_name]]))
    mean_per_row <- apply(rt_mat,1, row_mean)
    input <- lapply(col_names, function(y) as.data.frame(do.call(cbind, lapply(gc_peak_list_raw, function(x) x[y]))))
    input <- lapply(input, function(x){
        names(x) <- names(gc_peak_list_raw)
        x
    })

    input <- lapply(input, function(x){
        x <- cbind(mean_per_row, x)
        x
    })

    input <- lapply(input, function(x){
        names(x)[1] <- "mean_RT"
        x
    })
    names(input) <- col_names

    ### 7 Some documentation
    ### =====================
    # Call of align_chromatograms, List
    call <- as.list(match.call())[-1]
    # Defaults added to the function call List
    call <- function_call(call = call,FUN = align_chromatograms)
    Logbook[["Call"]] <- call
    ### =====================

    ### 8 Write output to text files
    ### ============================

    if (!is.null(write_output)) {
        # Take the text-file name as a prefix for each output file
        if (is.character(data)) {
            prefix <- strsplit(data,split = "/")
            prefix <- as.character(prefix[[1]][length(prefix[[1]])])
            prefix <- as.character(strsplit(prefix,split = ".txt"))
        } else {
            # For a List take "Aligned" as prefix
            prefix <- as.character(Logbook[["Call"]][["data"]])
        }
        file_names <- lapply(X = write_output, FUN = write_files,data = output, name = prefix)
        Logbook[["Output"]] <- file_names
    }

    ### 8 Documentation in Logbook
    ### ==========================

    Logbook[["Variation"]][["Input"]] <- unlist(align_var(gc_peak_list_raw,rt_col_name))
    Logbook[["Variation"]][["LinShift"]] <- unlist(align_var(gc_peak_list_linear,rt_col_name))
    Logbook[["Variation"]][["Aligned"]] <- unlist(align_var(gc_peak_list_aligned,rt_col_name))
    Logbook[["Date"]]["End"] <- as.character(strftime(Sys.time()))

    ### 9 Generate a list containing all returned output

    if ("GCalignR_Dummy" %in% col_names) {
        # remove the dummy variable when one was created
        gc_peak_list <- lapply(X = gc_peak_list, dummy_remove)
        gc_peak_list_aligned <- lapply(X = gc_peak_list_aligned, dummy_remove)
        gc_peak_list_linear <- lapply(X = gc_peak_list_linear, dummy_remove)
        gc_peak_list_raw <- lapply(X = gc_peak_list_raw, dummy_remove)
        output <- output[-2]
        input <- input[-2]
    }
    # substitute manipulated retention times with the input values
    output <- remove_linshifts(dx = output, rt_col_name = rt_col_name, Logbook = Logbook)
    gc_peak_list_aligned <- remove_linshifts2(dx = gc_peak_list_aligned, rt_col_name = rt_col_name, Logbook = Logbook)
    output_algorithm <- list(aligned = output,
                             heatmap_input = list(input_rts = rt_raw,
                                                  linear_transformed_rts = rt_linear,
                                                  aligned_rts = rt_aligned),
                             Logfile = Logbook,
                             aligned_list = gc_peak_list_aligned,
                             input_list = input_list,
                             input_matrix = input)

    class(output_algorithm) <- "GCalign"

    cat(paste0('\nAlignment completed!\n','Time:'),strftime(Sys.time(),format = "%Y-%m-%d %H:%M:%S"),'\n\n')
    return(output_algorithm)
}# end align_chromatograms
