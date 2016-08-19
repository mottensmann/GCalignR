#' Aligning Peaks based on retention times
#'
#'@description
#'\strong{align_chromatograms} is the core function of \code{\link{GCalignR}} that fasciliates the alignment of gas-chromatography data.
#' Read through the documentation here and take a look at the \href{../doc/GCalignR_step_by_step.html}{GCalignR Vignette}
#'
#' .
#'
#'@details
#' The Alignment of Peaks is archieved by running \strong{three major algorithms} always considering the complete
#' set of samples submitted to the function.
#' In brief: \strong{(1) Peaks are linearly shifted} to maximise similarity with a reference to account for systematic shifts in retention times caused by gas-chromatography processing. \strong{(2) Peaks of similar retention times are aligned} in order to match similar retention times to the same substance. During the algorithm proceeds, these clusters are continously revised and every peaks is moved to the optimal location(i.e. substance). \strong{(3) Peaks of similar retention time are merged} if they show smaller differences in mean retention times than expected by the achievable resolution of the gas-chromatography or the chemistry of the compounds are merged. This has to be specfied by the paramters \code{max_diff_peak2mean} and \code{min_diff_peak2peak}. Several optional processing steps are available, ranging from the removal of peaks representing contaminations (requires to include blanks as a control) to the removal of uninformative peaks that are present in just one sample.
#'
#'@param data
#' Two input formats are supported. The first option is the \strong{path to a text file} with extension \code{.txt} containing the gc-data. It is expected that the file is formatted following this principle: The first row contains sample names, the second row column names of the corresponding chromatograms.Starting with the third row, peak data are included, whereby matrices of single samples are concatenated horizontally. The matrix for each sample needs to consist of the same number of columns, at least two are required: The retention time and a measure of concentration (e.g. peak area or height). See the package vignette for an example. Alternatively the input may be a \strong{list of data frames}. Each data frame contains the peak data for a single individual with at least two variables, the retention time of the peak and the area under the peak. The variables need to have the same names across all samples (i.e. data frames). Also, each list element has to be named with the ID of the respective sample. See the vignette for an example. The data can be checked by running \code{\link{check_input}}.
#'
#'@param sep
#' The field separator character. The default is tab seperated (\code{sep = '\\t'}). See the "sep" argument in \code{\link[utils]{read.table}} for details.
#'
#'@param conc_col_name
#' Character naming a column used for the quantification of peaks (e.g. peak area or peak height).
#' The designated variable needs to be numeric.
#'
#'@param rt_col_name
#' Character string naming the column containing retention times.The designated variable needs to be numeric.
#'
#'@param write_output
#' Character vector of variables to write to a text file (e.g. \code{c("RT","Area")}.
#' The default is \code{NULL}. Names of the text files are concatenations of \code{data} and the output variables. Vector elements need to correspond to column names of \code{data}.
#'
#'@param rt_cutoff_low
#' Lower threshold under which retention times are cutted (i.e. 5 minutes). Default is NULL.
#'
#'@param rt_cutoff_high
#' Upper threshold above which retention times are cutted (i.e. 35 minutes). Default is NULL.
#'
#'@param reference
#' Character string of a sample to which all other samples are aligned to by means of a
#' linear shift (e.g. \code{"M3"}. The name has to correspond to an individual name given
#' in the first line of \code{data}. Alternatively a sample called \code{reference} can be included
#' in \code{data} containing user-defined peaks (e.g. an internal standard) to align the samples to. After the linear transformation the \code{reference} will be removed from the data.
#'
#'@param max_linear_shift
#' Defines a window to search for an optimal linear shift of chromatogram peaks with respect to the reference. Shifts are evaluated within - \code{max_linear_shift to + max_linear_shift}.The default is 0.05. The value of this parameter strongly depends on the linear trends the are expected for a given data set and should cover the expected range of linear shifts.
#'
#'@param max_diff_peak2mean
#' Defines the allowed deviation of retention times around the mean of the corresponding row (i.e. scored substance). Peaks with differing retention times are moved to a more appropriate mean retention time. The default is 0.03.
#'
#'@param min_diff_peak2peak
#' Defines the minimum difference in retention times among distinct substances. Substances that differ less, are merged if every sample contains either one or none of the respective compounds. This parameter is a major determinant in the classification of distinct peaks. Therefore careful consideration is required to adjust this setting to your needs (i.e. the resolution of your gas-chromatography pipeline). Large values may cause the merge of true peaks with similar retention times, if those are not simultaneously occuring within at least one individual. Especially for small sample sizes this could be the case just by chance. Small values can be considered as conservative and will reduce the number of scored substances. Default is 0.03.
#'@param blanks
#' Character vector of names of blanks. If specified, all substances found in any of the blanks will be removed from all samples, before the blanks are deleted from the aligned data. The names have to correspond to a name given in the first line of \code{data}.
#'
#'@param delete_single_peak
#' logical, determining whether substances that occur in just one sample are removed or not. By default single substances are retained in chromatograms.
#'
#'@param n_iter
#' integer indicating the number of iterations of the core alignment algorithm. Additional replications of the alignment and merging steps might be helpful to clean-up chromatograms, that otherwise show some remaining peak outliers mapped to the wrong mean retention time. Inspect alignment visually with a Heatmap \code{\link{gc_heatmap}}. Currently we recommend to run the algorithm once and check the output carefully to fine tune function paramters rather than running multiple iterations.
#'
#' @param merge_rare_peaks
#' logical determining whether peaks that are redundant by means of a similar retention times
#' are merged. If \code{TRUE} peaks are merged as long as only 5 % of the samples contain two peaks.
#' Always the peak with the higher abundance (i.e. peak area or peak height) is retained. Default is \code{FALSE}
#'
#'@return
#' Returns an object of class "GCalign" that is a a list with several objects that are listed below. Note, that the objects "heatmap_input" and "Logfile" are best inspected by calling the provided functions \emph{gc_heatmap} and \emph{print}.
#' \item{aligned}{Aligned gas-chromatography data subdivided into individual data frames for every variable. Samples are represented by columns, rows specifiy substances. The first column of every data frame is comprised of the mean retention time of the respective substance (i.e. row).}
#' \item{heatmap_input}{Data frames of retention times; used internally to create heatmaps}
#' \item{Logfile}{Included several lists summarizing the data; used to print diagonistics of the alignment}
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann (meinolf.ottensmann@@web.de)
#'
#'@import magrittr
#'
#'@examples
#' ## Load example data set
#' data("peak_data")
#' ## Subset for faster processing
#' peak_data <- peak_data[1:4]
#' peak_data <- lapply(peak_data,function(x) x <- x[1:80,])
#' out <- align_chromatograms(peak_data, conc_col_name = "area", rt_col_name = "time",
#'        rt_cutoff_low = 5, rt_cutoff_high = 45, reference = "M3",
#'          max_linear_shift = 0.05, max_diff_peak2mean = 0.03, min_diff_peak2peak = 0.03,
#'          blanks = NULL, delete_single_peak = TRUE)
#'@export
#'
align_chromatograms <- function(data, sep = "\t",conc_col_name=NULL, rt_col_name = NULL, write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.02, blanks = NULL,delete_single_peak = FALSE,n_iter=1,merge_rare_peaks=FALSE) {

## Check the format of the input
check_input(data,sep,write_output=write_output,blank=blanks)

## Preallocate a list to save a protocoll of the alignment
Logbook <- list()
Logbook[["Date"]]["Start"] <- as.character(strftime(Sys.time()))

## print start of GCalignR
cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%H:%M:%S")),'\n\n'))

## 1: Check function call
if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
if (is.null(conc_col_name)) stop("Column containing concentration of peaks is not specified. Define conc_col_name")
if (is.null(reference)) stop("Reference is missing. Specify a reference to align the others to")

## 2: Load Data
if (is.character(data)) { # Data is supplied as a Text File
    ind_names <- readr::read_lines(data, n_max = 1) %>% # Get Sample Names
        stringr::str_split(pattern = sep) %>%
        unlist()
    ind_names <- ind_names[ind_names != ""]
    col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>% # Get Variable Names
        stringr::str_split(pattern = sep) %>%
        unlist()
    col_names <- col_names[col_names != ""]
    col_names <- stringr::str_trim(col_names)
    ind_names <- stringr::str_trim(ind_names)

    gc_data <- utils::read.table(data, skip = 2, sep = sep, stringsAsFactors = F) # Get Peak Data
    gc_data <- gc_data[!(rowSums(is.na(gc_data)) == ncol(gc_data)), ] # Remove just NA-rows
    gc_data <- gc_data[,!(colSums(is.na(gc_data)) == nrow(gc_data))] # Remove empty rows
    gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))# Transform variables to numeric
    gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names) # convert to list

} else if (is.list(data)) { # Dat is supplied in a List
    col_names <- unlist(lapply(data, function(x) out <- names(x)))
    col_names <- names(data[[1]])
    ind_names <- names(data)
    gc_peak_list <- data
}

## Write some information about the input data to the Logfile
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

## 3: Start the processing of the gas-chromatography-data

## 3.1: Cut retention times
    gc_peak_list <- lapply(gc_peak_list, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_col_name)
    gc_peak_list_raw <- lapply(gc_peak_list, matrix_append, gc_peak_list)

## Write to the Logfile
if (!is.null(rt_cutoff_low) & is.null(rt_cutoff_high)){  # Lower threshold
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
}
if (!is.null(rt_cutoff_high) & is.null(rt_cutoff_low)){  # upper threshold
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
}
if (!is.null(rt_cutoff_high) & !is.null(rt_cutoff_low)){ # Lower and upper threshold
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- as.character(rt_cutoff_low)
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- as.character(rt_cutoff_high)
}
if(is.null(rt_cutoff_high) & is.null(rt_cutoff_low)){#None
    Logbook[["Filtering"]]["RT_Cutoff_Low"] <- "None"
    Logbook[["Filtering"]]["RT_Cutoff_High"] <- "None"
}

## 3.2: Linear Transformation of peak retention times

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

## 3.3 Align peaks & Merge Rows

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

## 4. Clean Up Chromatograms
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

## 5: Create Output for Heatmaps
    rt_raw <- rt_extract(gc_peak_list = gc_peak_list_raw,rt_col_name =rt_col_name) # Initial input
    rt_linear <- rt_extract(gc_peak_list = gc_peak_list_linear,rt_col_name =rt_col_name) # after linear shifts
    rt_aligned <- rt_extract(gc_peak_list = gc_peak_list_aligned,rt_col_name =rt_col_name) # final output

    row_mean <- function(x) {if (all(x==0)) 0 else mean(x[x!=0])} # mean per row without 0
    mean_per_row <- apply(rt_mat,1, row_mean)

## 6: Create output matrices for Variables
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

## 7: Write output to text files

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

# 8: Documentation in the Logfile

    Logbook[["Variation"]][["Input"]] <- unlist(align_var(gc_peak_list_raw,rt_col_name))
    Logbook[["Variation"]][["LinShift"]] <- unlist(align_var(gc_peak_list_linear,rt_col_name))
    Logbook[["Variation"]][["Aligned"]] <- unlist(align_var(gc_peak_list_aligned,rt_col_name))
    Logbook[["Date"]]["End"] <- as.character(strftime(Sys.time()))

    call <- as.list(match.call())[-1] # Call of align_chromatograms, List
    call <- function_call(call = call,FUN = align_chromatograms) # Defaults added
    Logbook[["Call"]] <- call

## 9 Generate output of the function call as a list

        output_algorithm <- list(aligned=output, #summary=align_summary,
                             heatmap_input=list(input_rts=rt_raw,
                             linear_transformed_rts=rt_linear,aligned_rts=rt_aligned),
                             Logfile=Logbook)

    class(output_algorithm) <- "GCalign" # Define GCalignR object calss

    cat(paste('Alignment was Successful!\n','Time:'),strftime(Sys.time(),format = "%H:%M:%S"),'\n')
    return(output_algorithm)
}
