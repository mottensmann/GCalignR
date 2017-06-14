#' Aligning peaks based on retention times
#'
#'@description
#' This is the core function of \code{\link{GCalignR}} to align peak data. The input data is a peak list.
#' Read through the documentation below and take a look at the vignettes for a thorough introduction.
#'
#'@details
#' The alignment of peaks is achieved by running \strong{three major algorithms} always considering
#' the complete set of samples submitted to the function.
#' In brief: \strong{(1) Chromatograms (more correctly, their peaks) are shifted} to maximise
#' similarity with a reference to account for systematic shifts in retention times
#' caused by gas-chromatography processing. \strong{(2) Peaks of similar retention times are aligned}
#' in order to match similar retention times to the same substance. During the algorithm proceeds,
#' these clusters are continuously revised and every peaks is moved to the optimal
#' location(i.e. substance). \strong{(3) Peaks of similar retention time are merged} if
#' they show smaller differences in mean retention times than expected by the achievable
#' resolution of the gas-chromatography or the chemistry of the compounds are merged.
#' This has to be specfied by the parameters \code{max_diff_peak2mean} and \code{min_diff_peak2peak}.
#' Several optional processing steps are available, ranging from the removal of peaks representing
#' contaminations (requires to include blanks as a control) to the removal of uninformative peaks
#' that are present in just one sample.
#'
#'@param data
#' Chemical data that has to be aligned. All variables that are available need to be included in order to align these measures based on the retention time, which is the only mandatory variable.
#' Two input formats are supported. The first option is the \strong{path to a plain text file} with extension ".txt" containing the gc-data. It is expected that the file is formatted following this
#' principle: The first row contains sample names, the second row column names of the corresponding
#' chromatograms. Starting with the third row, peak data are included, whereby matrices of single
#' samples are concatenated horizontally (see the vignette or example data). The matrix for each
#' sample needs to consist of the same number of columns, at least one is required that contains the retention times of peaks. See the \href{../doc/GCalignR_step_by_step.html}{vignette} for an example. Alternatively the input may be a \strong{list of data frames}. Each data frame contains
#' the peak data for a single individual with at least one column of containing retention times of peaks. Variables need to have the same names across all samples (i.e. data frames). Also, each list element has to be named with the ID of the respective sample.
#' The format can be checked by running \code{\link{check_input}}.
#'
#'@param sep
#' The field separator character. The default is tab separated (\code{sep = '\\t'}).
#' See the "sep" argument in \code{\link[utils]{read.table}} for details.
#'
#'@param rt_col_name
#' Character string - the name of the column containing the retention times.The variable needs to
#' be numeric and the decimal separator needs to be a point.
#'
#'@param write_output
#' Character vector of variables to write to a text file (e.g. \code{c("RT","Area")}.
#' Vector elements need to correspond to column names of \code{data}. Writing output is optional.
#'
#'@param rt_cutoff_low
#' Lower threshold under which retention times are removed (i.e. 5 minutes). Default NULL.
#'
#'@param rt_cutoff_high
#' Upper threshold above which retention times are removed (i.e. 35 minutes). Default NULL.
#'
#'@param reference
#' Character string of a sample to which all other samples are aligned by means of a
#' linear shift (e.g. \code{"M3"}. The name has to correspond to an individual name given
#' in the first line of \code{data}. Alternatively a sample called \code{reference} can be included
#' in \code{data} containing user-defined peaks (e.g. an internal standard) to align the samples to.
#' After the linear transformation the \code{reference} will be removed from the data.
#'
#'@param max_linear_shift
#' This value defines the maximum time that one chromatogram is expected to be deviating in retention times
#' from another chromatogram. To correct for these systematic shifts, the algorithm potentially adds the same
#' retention time to all peaks within a chromatogram to maximise the number of shared peaks with
#' the reference. We recommend to start with the default of 0.02 (minutes) and increase if necessary.
#'
#'@param max_diff_peak2mean
#' Numeric value defining the allowed deviation of the retention time of a given peak from the mean of the corresponding row (i.e. scored substance). Defaults to 0.02 (minutes). This parameter reflects the retention time range in which peaks across samples are still counted as the same 'substance',
#' i.e. sorted in one row.
#'
#'@param min_diff_peak2peak
#' Numeric values defining the expected minimum difference in retention times among different substances.
#' Retention time rows that differ less, are therefore merged if every sample contains either
#' one or none of the respective compounds.
#' This parameter is a major determinant in the classification of distinct peaks.
#' Therefore careful consideration is required to adjust this setting to your needs
#' (e.g. the resolution of your gas-chromatography pipeline).
#' Large values may cause the merge truly different substances with similar retention times, if those are not
#' simultaneously occurring within at least one individual, which might occur by chance for small sample
#' sizes. It is therefore recommended to set the value much lower (e.g. 0.02) when few individuals are
#' analysed. Defaults to 0.08 (minutes).
#'
#'@param blanks
#' Character vector of names of negative controls. Substances found in any of the blanks will be
#' removed from all samples, before the blanks are deleted from the aligned data.
#'
#'@param delete_single_peak
#' Logical, determining whether substances that occur in just one sample are removed or not. Default FALSE.
#'
#'@return
#' Returns an object of class "GCalign" that is a a list containing several objects that are listed below.
#' Note, that the objects "heatmap_input" and "Logfile" are best inspected by calling the provided functions \emph{gc_heatmap} and \emph{print}.
#' \item{aligned}{Aligned gas-chromatography data subdivided into individual data frames for every variable.
#' Samples are represented by columns, rows specify substances. The first column of every data frame
#' is comprised of the mean retention time of the respective substance (i.e. row). The aligned data
#' can be used for further statistical analyses in other packages and also directly written
#' to .txt file by specifying the write_output argument in align_chromatograms}
#' \item{heatmap_input}{Data frames of retention times; used internally to create heatmaps}
#' \item{Logfile}{Includes several lists summarizing the data; used to print diagnostics of the alignment}
#' \item{input_list}{List of data frames. Data frames are comprised of the raw of a sample prior to aligning}
#' \item{input_matrix}{List of data frames. Data frames are matrices of input variables}
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
#' ## align data
#' out <- align_chromatograms(peak_data, rt_col_name = "time",
#' rt_cutoff_low = 10, rt_cutoff_high = 30, reference = "M2",
#' max_linear_shift = 0.02)
#'
#'@export
#'
align_chromatograms <- function(data, sep = "\t", rt_col_name = NULL,
    write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
    max_linear_shift = 0.02, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.08, blanks = NULL,
    delete_single_peak = FALSE) {

# Print start
cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%Y-%m-%d %H:%M:%S")),'\n\n'))

# Iteration have been deleted as function paramter.
iterations = 1

### 1. Checks preparations
### ======================

# 1.1 Stop execution if mandatory checks are not passed
if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
x <- check_input(data,sep,write_output = write_output,blank = blanks,reference = reference,rt_col_name = rt_col_name)
if (x != TRUE) stop("Processing not possible: check warnings below and change accordingly in order to proceed")


# 1.2 Create a "Logbook" to record alignment steps and parameters
Logbook <- list()
Logbook[["Date"]]["Start"] <- as.character(strftime(Sys.time()))

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
    gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F)
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

} else if (is.list(data)) { # data is in a list
    col_names <- unlist(lapply(data, function(x) out <- names(x)))
    col_names <- names(data[[1]])
    ind_names <- names(data)
    gc_peak_list <- lapply(data,matrix_append,gc_peak_list = data, val = "NA") # same dimensions of dfs

} # end load data

# save gc_peak_list for documentation purposes
input_list <- gc_peak_list

# Write some information about the input data to the Logfile
if (!is.null(reference)) {
    if (reference == "reference") {
    cat(paste0('Data for ',as.character(length(ind_names) - 1),' samples loaded.'))
    Logbook[["Input"]]["Samples"] <- length(ind_names) - 1
    }
    } else {
    cat(paste0('Data for ',as.character(length(ind_names)),' samples loaded.'))
    Logbook[["Input"]]["Samples"] <- length(ind_names)
}
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
if (reference == "reference") {
    # reference is not a true sample and is used exclusevely for linear alignment
    cat(paste0('\nStart correcting linear shifts with ',"\"",as.character(reference),"\"",' as a reference ...'))

gc_peak_list_linear <- linear_transformation(gc_peak_list = gc_peak_list, max_linear_shift = max_linear_shift, step_size = 0.01, reference = reference, rt_col_name = rt_col_name, Logbook = Logbook)
    Logbook <- gc_peak_list_linear[["Logbook"]]
    gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
    gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
    gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames,col_names)
    # remove reference after the systematic errors were corrected
    gc_peak_list_linear <- gc_peak_list_linear[-which(names(gc_peak_list_linear) == reference)]
    # remove the reference from the input retention time matrix
    gc_peak_list_raw <- gc_peak_list_raw[-which(names(gc_peak_list_raw) == reference)]
}  else {
    gc_peak_list_linear <- linear_transformation(gc_peak_list = gc_peak_list, max_linear_shift = max_linear_shift, step_size = 0.01, reference = reference, rt_col_name = rt_col_name, Logbook = Logbook)
    Logbook <- gc_peak_list_linear[["Logbook"]]
    gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
    gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
    gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames, col_names)
}
} else if (is.null(reference)) {
    # If no reference was specified by the user, the reference is determined, such
    # that the sample with the highest avarage similarity to all other samples is used.
    cat("\nNo reference was specified. Hence, a reference is selected automatically ... ")
    best.ref <- choose_optimal_reference(gc_peak_list = gc_peak_list, rt_col_name = rt_col_name)
    # set the reference
    reference <- best.ref[["sample"]]
    cat("done\n")
    text <- paste0("'",reference,"' was selected on the basis of highest average similarity to all samples (score = ",best.ref[["score"]],").")
    cat(stringr::str_wrap(paste(text,collapse = ""),width = 80,exdent = 0,indent = 0))
}# new end
     cat(paste0('\nStart correcting linear shifts with ',"\"",as.character(reference),"\"",' as a reference ...'))
    gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift = max_linear_shift,
        step_size = 0.01, reference = reference, rt_col_name = rt_col_name, Logbook = Logbook)
    Logbook <- gc_peak_list_linear[["Logbook"]]
    gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
    gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
    gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames,col_names)
#} old end

    Logbook[["Input"]]["Reference"] <- reference
    # why is there a cat?
    cat(" done\n")
    # equalise chromatograms sizes
    gc_peak_list_linear <- lapply(gc_peak_list_linear, matrix_append, gc_peak_list_linear)

# 3.3 Align peaks

    cat(c('Start aligning peaks ... ','this might take a while!\n'))

    # create corresponding lists
    gc_peak_list_aligned <- gc_peak_list_linear
    no_peaks <- matrix(NA,nrow = iterations,ncol = 1)
    merged_peaks <- matrix(NA, nrow = iterations,ncol = 1)

    # Iterations are currently limited to one loop
for (R in 1:iterations) {
    gc_peak_list_aligned <- align_peaks(gc_peak_list_aligned, max_diff_peak2mean = max_diff_peak2mean,
        iterations = iterations, rt_col_name = rt_col_name,R = R)
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
    cat("\nMerge redundant substances ... ")
    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned,
        min_diff_peak2peak = min_diff_peak2peak, rt_col_name = rt_col_name)
    # estimate Number of merged peaks
    merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])
    # Number of peaks after merging
    no_peaks[R] <- nrow(gc_peak_list_aligned[[1]])

    Logbook[["Aligned"]]["total"] <- no_peaks[R]

    cat('done')
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
    # delete empty rows again
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)
}#end iterative loop of aligning & merging

cat(paste0('\n','Alignment completed'),'\n\n')

### 4 Cleaning chromatograms
### ========================
    gc_peak_list_aligned <- gc_peak_list_aligned[match(names(gc_peak_list_raw),names(gc_peak_list_aligned))]

    # delete peaks present in blanks, then remove the blank itself
if (!is.null(blanks)) {
    cat("Remove contaminations and remove blanks ... ")
    # delete peaks present in blanks
    delete_blank <- function(blank, gc_peak_list_aligned) {
        # indices of peaks
    del_substances <- which(gc_peak_list_aligned[[blank]][[rt_col_name]] > 0)
    # remove peaks, delete row
    chroma_out <- lapply(gc_peak_list_aligned, function(x) x[-del_substances,])
    # removes the blank from the list
    chroma_out[blank] <- NULL
    return(chroma_out)
    }
    # Number of Peaks including Blanks
    N <- nrow(gc_peak_list_aligned[[1]])
    # delete all blanks
for (i in blanks) {gc_peak_list_aligned <- delete_blank(i, gc_peak_list_aligned)}
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

### 7 Some protocollation
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
    output_algorithm <- list(aligned = output,
                             heatmap_input = list(input_rts = rt_raw,
                             linear_transformed_rts = rt_linear,aligned_rts = rt_aligned),
                             Logfile = Logbook, aligned_list = gc_peak_list_aligned ,input_list = input_list, input_matrix = input)

    class(output_algorithm) <- "GCalign"

    cat(paste0('\nAlignment was successful!\n','Time:'),strftime(Sys.time(),format = "%Y-%m-%d %H:%M:%S"),'\n\n')
    return(output_algorithm)
}# end align_chromatograms
