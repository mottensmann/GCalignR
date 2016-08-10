#' Aligning Peaks based on retention times
#'
#'@description
#'\code{align_chromatograms()} is the core function of \code{\link{GCalignR}}. Check out the vignette
#'to get started: utils::browseVignettes(package = "GCalignR")
#'
#'@details
#'The Alignment of Peaks is archieved by running three major algorithms always considering the complete
#'set of samples submitted to the function. In brief: All peaks are linearly shifted to
#'maximise similarity with a reference to account for systematic shifts in retention times
#'caused by gas-chromatography processing. Second, peaks of similar retention times are
#'transfered step-wise to the same location (i.e row) in order to group substances. In a third step peaks (i.e.substances)
#'that show smaller differences in mean retention times than expected by the achievable resolution
#'of the gas-chromatography or the chemistry of the compounds are merged. Several optional processing steps are available
#'ranging from the removal of peaks representing contaminations (requires to include blanks as a control) to the removal
#'of uninformative peaks that are present in just one sample. Further Details can be found in
#'
#'@param data
#' Two options. (1) Path to a file with extension \code{.txt} containing the gc-data. It is expected that the
#' file is formatted following this principle: The first row contains sample names, the second row column names of the corresponding chromatograms.
#' Starting with the third row, peak data are included, whereby matrices of single samples are concatenated horizontally. The
#' matrix for each sample needs to consist of the same number of columns, at least two are required:The
#' retention time and a measure of concentration (e.g. peak area or height). See the package
#' vignette for an example. (2) The alternative input is a list of data.frames. Each data.frame
#' contains the peak data for a single individual with at least two variables, the retention time of
#' the peak and the area under the peak. The variables have to have the same name across all samples
#' (data.frames). Also, each list element (i.e. each data.frame) has to be named with the ID of
#' the individual. See the vignette for an example. The data can be checked by running \code{\link{check_input}}
#'
#'@param sep
#'The field separator character. Values on each line of the file are separated by this
#'character. The default is tab seperated (\code{sep = '\\t'}). See the \code{sep} argument in \code{\link[utils]{read.table}} for details.
#'
#'@param conc_col_name
#'Character string naming a column used for the quantification of peaks (e.g. peak area or peak height).
#'The designated variable needs to be numeric.
#'
#'@param rt_col_name
#'Character string naming the column containing retention times.The designated variable needs to be numeric.
#'
#'@param write_output
#' Character vector of variables to write to a text file (e.g. \code{c("RT","Area")}.
#' The default is \code{NULL}.Names of the text files are concatenations of \code{datafile} and the output variables.
#' Strings need to correspond to column names of \code{data}.
#'
#'@param rt_cutoff_low
#'Lower threshold under which retention times are cutted (i.e. 5 minutes). Default is NULL.
#'
#'@param rt_cutoff_high
#'Upper threshold above which retention times are cutted (i.e. 35 minutes). Default is NULL.
#'
#'@param reference
#'Character string of a sample to which all other samples are aligned to by means of a
#' linear shift (e.g. \code{"individual3"}. The name has to correspond to an individual name given
#' in the first line of \code{data}. Alternatively a sample called \code{reference} can be included
#' in \code{data} containing user-defined peaks (e.g. an internal standard) to align the samples to. After the linear
#' transformation \code{reference} will be removed from the data.
#'
#'@param max_linear_shift
#' Defines a window to search for an optimal linear shift of chromatogram peaks
#' with respect to the reference. Shifts are evaluated within - \code{max_linear_shift to + max_linear_shift}.
#' Default is 0.05.
#'
#'@param max_diff_peak2mean
#' Defines the allowed deviation of retention times around the mean of
#' the corresponding row (i.e. substance). Peaks with differing retention times are moved
#' to a more appropriate mean retention time. Default is 0.02.
#'
#'@param min_diff_peak2peak
#'Defines the minimum difference in retention times among distinct
#'substances. Substances that differ less, are merged if every sample contains either one
#'or none of the respective compounds. This parameter is a major determinant in the classification
#'of distinct peaks. Therefore careful consideration is required to adjust this setting to
#'your needs (i.e. the resolution of your gas-chromatography pipeline). Large values may cause the
#'merge of true peaks with similar retention times, if those are not simultaneously occuring within at least
#'one individual. Especially for small sample sizes this could be the case just by chance
#'Small values can be considered as conservative. Default is 0.02.
#'
#'@param blanks
#'Character vector of names of blanks. If specified, all substances found in any of the blanks
#'will be removed from all samples, before the blanks are deleted from the aligned data.
#'The names have to correspond to a name given in the first line of \code{data}.
#'
#'@param delete_single_peak
#'logical, determining whether substances that occur in just one sample are
#'removed or not. By default single substances are retained in chromatograms.
#'
#'@param n_iter
#' integer indicating the number of iterations of the core alignment algorithm. Additional replications of
#' the alignment and merging steps might be helpful to clean-up chromatograms, that otherwise show
#' some remaining peak outliers mapped to the wrong mean retention time. Inspect alignment visually
#' with a Heatmap \code{\link{gc_heatmap}}.
#'
#' @param merge_rare_peaks
#' logical determining whether peaks that are redundant by means of a similar retention times
#' are merged. If \code{TRUE} peaks are merged as long as only 5 % of the samples contain two peaks.
#' Always the peak with the higher abundance (i.e. peak area or peak height) is retained.
#'
#'
#'@return
#' Returns an object of class GCalign that is a a list with the following elements:
#' \item{call}{function call}
#' \item{chroma_aligned}{a list containing data.frames with the aligned variables}
#' \item{rt_raw}{a data.frame with the retention times before alignment}
#' \item{rt_linear}{a data.frame with the retention times after the linear transformation}
#' \item{rt_aligned}{a data.frame with the final aligned retention times}
#' \item{align_summary}{a data.frame}{A basic summary of the alignment process}
#'
#'  @author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#' @import magrittr
#'
#' @examples
#' data(seal_peaks)
#' seal_peaks <- seal_peaks[1:4]
#' out <- align_chromatograms(seal_peaks, conc_col_name = "area", rt_col_name = "time",
#'        rt_cutoff_low = 5, rt_cutoff_high = 45, reference = "M3",
#'          max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03,
#'          blanks = NULL, delete_single_peak = TRUE)
#'
#' @export
#'
align_chromatograms <- function(data, sep = "\t",conc_col_name=NULL, rt_col_name = NULL, write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.02, blanks = NULL,delete_single_peak = FALSE,n_iter=1,merge_rare_peaks=FALSE) {

check_input(data,sep,write_output=write_output,blank=blanks)

Logbook <- list() # List saving the main workflow
Logbook[["Date"]]["Start"] <- as.character(strftime(Sys.time()))

cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%H:%M:%S")),'\n\n')) # print time of GCaligner Start

# 1: Check function call
# ----------------------
if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
if (is.null(conc_col_name)) stop("Column containing concentration of peaks is not specified. Define conc_col_name")
if (is.null(reference)) stop("Reference is missing. Specify a reference to align the others to")

# 2: Load Data
#-------------
if (is.character(data)) {# TextFile
    ind_names <- readr::read_lines(data, n_max = 1) %>% # Sample Names
        stringr::str_split(pattern = sep) %>%
        unlist()
    ind_names <- ind_names[ind_names != ""]

    col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>% # Variable Names
        stringr::str_split(pattern = sep) %>%
        unlist()
    col_names <- col_names[col_names != ""]
    col_names <- stringr::str_trim(col_names)
    ind_names <- stringr::str_trim(ind_names)

    gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F) # Peak Data
    gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ] # Remove just NA-rows
    gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))] # Remove empty rows
    gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))# Transform variables to numeric
    gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names) # convert to list

} else if (is.list(data)) {#List
    col_names <- unlist(lapply(data, function(x) out <- names(x)))
    col_names <- names(data[[1]])
    ind_names <- names(data)
    gc_peak_list <- data
}#End LoadData

if(reference == "reference"){
    cat(paste0('GC-data for ',as.character(length(ind_names) - 1),' samples loaded\n'))
    Logbook[["Input"]]["Samples"] <- length(ind_names) - 1
} else {
    cat(paste0('GC-data for ',as.character(length(ind_names)),' samples loaded\n'))
    Logbook[["Input"]]["Samples"] <- length(ind_names)
}
    Logbook[["Input"]]["Range"] <- paste((range(peak_lister(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name))),collapse = "-")
    Logbook[["Input"]]["File"] <- as.character(as.character(match.call()["data"]))
    Logbook[["Input"]]["Reference"] <- reference
    Logbook[["Input"]]["Retention_Time"] <- rt_col_name
    Logbook[["Input"]]["Concentration"] <- conc_col_name
    Logbook[["Input"]]["Peaks"] <- peak_counter(gc_peak_list = gc_peak_list,rt_col_name = rt_col_name)
if(!is.null(blanks)) Logbook[["Input"]][["Blanks"]] <- paste(blanks,collapse = "; ") # Only created if blanks!=NULL

# 3: Processing
#---------------

# 3.1: Cut retention times
#-------------------------
    gc_peak_list <- lapply(gc_peak_list, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_col_name)
    gc_peak_list_raw <- lapply(gc_peak_list, matrix_append, gc_peak_list)

if (!is.null(rt_cutoff_low) & is.null(rt_cutoff_high)){ #Low
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
}
if (!is.null(rt_cutoff_high) & is.null(rt_cutoff_low)){#High
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
}
if (!is.null(rt_cutoff_high) & !is.null(rt_cutoff_low)){#Low+High
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
}
if(is.null(rt_cutoff_high) & is.null(rt_cutoff_low)){#None
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
}

# 3.2: Linear Transformation
#---------------------------

    # Round retention times to two decimals
    round_rt <- function(gc_peak_df,rt_col_name=rt_col_name){
    gc_peak_df[rt_col_name] <- round(gc_peak_df[rt_col_name],digits = 2)
    return(gc_peak_df)
    }

    gc_peak_list <- lapply(X = gc_peak_list,FUN = round_rt)

    cat(paste0('\nStart Linear Transformation with ',"\"",as.character(reference),"\"",' as a reference ...'))

if(reference=="reference"){ # reference is not a true sample and is used exclusevely for linear alignment
    gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift=max_linear_shift, step_size=0.01,
                                                     reference = reference, rt_col_name = rt_col_name,
                                                     Logbook = Logbook)
    Logbook <- gc_peak_list_linear[["Logbook"]]
    gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
    gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
    gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames,col_names)
    gc_peak_list_linear <- gc_peak_list_linear[-which(names(gc_peak_list_linear)==reference)] # remove reference
    gc_peak_list_raw <- gc_peak_list_raw[-which(names(gc_peak_list_raw)==reference)] # remove the reference

}else{
    gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift=max_linear_shift, step_size=0.01,
                                                     reference = reference, rt_col_name = rt_col_name,
                                                     Logbook = Logbook)
    Logbook <- gc_peak_list_linear[["Logbook"]]
    gc_peak_list_linear <- gc_peak_list_linear[["chroma_aligned"]]
    gc_peak_list_linear <- lapply(gc_peak_list_linear,function(x) data.frame(x))
    gc_peak_list_linear <- lapply(gc_peak_list_linear, correct_colnames,col_names)
}
    cat("Done\n\n")

    gc_peak_list_linear <- lapply(gc_peak_list_linear, matrix_append, gc_peak_list_linear) # equalise chromatograms sizes

# 3.3 Align peaks & Merge Rows
#-----------------------------

    cat(c('Start Alignment of Peaks ... ','This might take a while!\n','\n'))

    Fun_Fact()# Random trivia
    gc_peak_list_aligned <- gc_peak_list_linear # avoid to overwrite the list of linear transformed peaks
    no_peaks <- matrix(NA,nrow = n_iter,ncol = 1) #
    merged_peaks <- matrix(NA, nrow = n_iter,ncol = 1)

for (R in 1:n_iter){ # Allows to iteratively execute the algorithm
    gc_peak_list_aligned <- align_peaks(gc_peak_list_aligned, max_diff_peak2mean = max_diff_peak2mean, n_iter = n_iter, rt_col_name = rt_col_name,R=R)
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)# mean rt per row
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) { # remove empty rows
    keep_rows <- which(!is.na(average_rts))
    out <- x[keep_rows, ]
    }) # Alignment is done

    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name) # mean rt per row
    no_peaks[R] <- nrow(gc_peak_list_aligned[[1]]) # N before merging
    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name)# merge rows
    merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])#estimate N of merged peaks
    no_peaks[R] <- nrow(gc_peak_list_aligned[[1]]) # N after merging
    Logbook[["Aligned"]]["Peaks"] <- no_peaks[R] # After the merging

if(R==n_iter & merge_rare_peaks==TRUE){ # allows to merge rows where a few individuals carry two peaks (<5%)
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name,criterion="proportional",conc_col_name = conc_col_name)
    merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])
    rare_peak_pairs <- 0 # cause they were eliminated
}

    cat('\n\nMerged Redundant Peaks')
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)     # delete empty rows again

} # End of iterative alignment and merging

    cat(paste0('\n\n','Peak Alignment Done'),'\n\n')

# 4. Clean Up Chromatograms
#--------------------------
    gc_peak_list_aligned <- gc_peak_list_aligned[match(names(gc_peak_list_raw),names(gc_peak_list_aligned))]

if (!is.null(blanks)) { # delete peaks present in blanks, then remove the blanks
    delete_blank <- function(blank, gc_peak_list_aligned) { # delete peaks present in blanks
    del_substances <- which(gc_peak_list_aligned[[blank]][[rt_col_name]] > 0) # indices of peaks
    chroma_out <- lapply(gc_peak_list_aligned, function(x) x[-del_substances,])# remove peaks, delete row
    chroma_out[blank] <- NULL # removes the blank from the list
    return(chroma_out)
}
    N <- nrow(gc_peak_list_aligned[[1]]) # Number of Peaks including Blanks
for (i in blanks) {gc_peak_list_aligned <- delete_blank(i, gc_peak_list_aligned)} # delete all blanks
    N <- N -  nrow(gc_peak_list_aligned[[1]])
    Logbook[["Filtering"]]["Blank_Peaks"] <- N
    Logbook[["Aligned"]]["In_Blanks"] <- N
    cat('Blank Peaks deleted & Blanks removed\n\n')
}

    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]])) # RT2Matrix
    singular_peaks <- nrow(rt_mat) # To estimate the number of deleted peaks
if (delete_single_peak) { # find single retention times in rows
    single_subs_ind <- which(rowSums(rt_mat > 0) == 1)
    if (length(single_subs_ind) > 0) {             # delete substances occuring in just one individual
            gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) x[-single_subs_ind, ])
    }
    cat(paste('Single Peaks deleted:',as.character(length(single_subs_ind)),'have been removed\n\n'))
    Logbook[["Aligned"]]["Singular"] <- length(single_subs_ind)
}
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]])) # final RT2Matrix
    singular_peaks <- singular_peaks - nrow(rt_mat) # how many were deleted, if any
    Logbook[["Aligned"]]["Retained"] <- nrow(gc_peak_list_aligned[[1]])

# 5: Create Output for Heatmaps
#------------------------------
    rt_raw <- rt_extract(gc_peak_list = gc_peak_list_raw,rt_col_name =rt_col_name) # Initial input
    rt_linear <- rt_extract(gc_peak_list = gc_peak_list_linear,rt_col_name =rt_col_name) # after linear shifts
    rt_aligned <- rt_extract(gc_peak_list = gc_peak_list_aligned,rt_col_name =rt_col_name) # final output

    row_mean <- function(x) {if (all(x==0)) 0 else mean(x[x!=0])} # mean per row without 0
    mean_per_row <- apply(rt_mat,1, row_mean)

# 6: Create output matrices for Variables
#----------------------------------------
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

# 7: Write to txt files
#----------------------

if (!is.null(write_output)){
    if(is.character(data)){ # Take the text-file name as a prefix for each output file
        prefix <- strsplit(data,split = "/")
        prefix <- as.character(prefix[[1]][length(prefix[[1]])])
        prefix <- as.character(strsplit(prefix,split = ".txt"))
    } else {
        prefix <- "Aligned" # For a List take "Aligned" as prefix
    }
        write_files <- function(x) {
        utils::write.table(output[[x]],
               file = paste0(prefix,"_",x, ".txt"), sep = "\t", row.names = FALSE)
    }
        Logbook[["Output"]] <- lapply(write_output,function(x) paste0(prefix,"_",x, ".txt"))
        lapply(write_output, write_files)
}

# 8: Documentation
#-----------------
    Logbook[["Variation"]][["Input"]] <- unlist(align_var(gc_peak_list_raw,rt_col_name))
    Logbook[["Variation"]][["LinShift"]] <- unlist(align_var(gc_peak_list_linear,rt_col_name))
    Logbook[["Variation"]][["Aligned"]] <- unlist(align_var(gc_peak_list_aligned,rt_col_name))
    Logbook[["Date"]]["End"] <- as.character(strftime(Sys.time()))

    call <- as.list(match.call())[-1] # Call of align_chromatograms
    call <- function_call(call = call,FUN = align_chromatograms) # Defaults added
    Logbook[["Call"]] <- call
    output_algorithm <- list(aligned=output, #summary=align_summary,
                             heatmap_input=list(input_rts=rt_raw,
                             linear_transformed_rts=rt_linear,aligned_rts=rt_aligned),
                             Logfile=Logbook)

    class(output_algorithm) <- "GCalign" # name of list

    cat(paste('Alignment was Successful!\n','Time:'),strftime(Sys.time(),format = "%H:%M:%S"),'\n')
    return(output_algorithm)
}#End Function
