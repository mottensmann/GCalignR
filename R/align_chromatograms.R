#' Aligning chromatograms based on retention times
#'
#'@description
#'\code{align_chromatograms()} is the core function of \code{\link{GCalignR}}.
#'
#'@details
#'Alignment is archieved by running three major algorithms always considering the complete
#'set of samples submitted to the function. In brief: All chromatograms are linearly shifted to
#'maximise similarity with a common reference to account for systematic shifts in retention times
#'caused by gas-chromatography processing. Second, peaks of similar retention times are step-wise
#'transfered to the same location (i.e row) in order to group substances. In a third step peaks (i.e.substances)
#'that show smaller differences in mean retention times than expected by the achievable resolution
#'of the gas-chromatography or the chemistry of the compounds are merged. Optional processing includes the
#'removal of peaks present in blanks (i.e. contaminations) and peaks that a uniquely found within just
#'one sample.
#'
#'@param data
#' Two options. (1) The most common input is the path of a file with extension \code{.txt} to load data from.
#' The first row needs to contain sample names, the second row column names of the corresponding chromatograms. Starting with the
#' third row, peak data are included, whereby matrices of single samples are concatenated horizontally. The
#' matrix for each sample needs to consist of the same number of columns, at least two are required:The
#' retention time and a measure of concentration (e.g. peak area or height). See the package
#' vignette for an example. (2) The alternative input is a list of data.frames. Each data.frame
#' contains the peak data for a single individual with at least two variables, the retention time of
#' the peak and the area under the peak. The variables have to have the name name across all samples
#' (data.frames). Also, each list element (i.e. each data.frame) has to be named with the ID of
#' the individual. See the vignette for an example.
#'
#'@param sep
#'The field separator character. Values on each line of the file are separated by this
#'character. The default is tab seperated (sep = '\\t'). See \code{sep} argument in \code{\link[utils]{read.table}} for details.
#'
#'@param conc_col_name
#'Character string naming a column used for quantification of peaks (e.g. peak area or peak height).
#'The designated variable needs to be numeric.
#'
#'@param rt_col_name
#'Character string naming the column holding retention times.The designated variable needs to be numeric.
#'
#'@param write_output
#' Character vector of variables to write to a text file (e.g. \code{c("RT","Area")}.
#' The default is \code{NULL}.Names of the text files are concatenations of \code{datafile} and the output variables.
#' Strings need to correspond to column name of the \code{datafile}.
#'
#'@param rt_cutoff_low
#'Lower threshold under which retention times are cutted (i.e. 5 minutes). Default is NULL.
#'
#'@param rt_cutoff_high
#'Upper threshold above which retention times are cutted (i.e. 35 minutes). Default is NULL.
#'
#'@param reference
#'Character string of a sample to which all other samples are aligned to by means of a
#'  linear shift (e.g. "individual3"). The name has to correspond to an individual name given
#'  in the first line of \code{datafile}.
#'
#'@param max_linear_shift
#' Defines a window to search for an optimal linear shift of chromatogram peaks
#' with respect to the reference. Shifts are evaluated within - \code{max_linear_shift to + max_linear_shift}.
#' Default is 0.05.
#'
#'@param max_diff_peak2mean
#' Defines the allowed deviation of retention times around the mean of
#' the corresponding row. Peaks with differing retention times are moved
#' to the appropriate mean retention time. Default is 0.02.
#'
#'@param min_diff_peak2peak
#'Defines the desired minimum difference in retention times among distinct
#'substances. Substances that differ less, are merged if every sample contains either one
#'or none of the respective compounds. This parameter is a major determinant in the classification
#'of distinct peaks. Therefore careful consideration is required to adjust this setting to
#'your needs (i.e. the resolution of your gas-chromatography pipeline). Large values might cause the
#'merge of true peaks, if those are not occuring within one individual, which might happen by chance
#'with a low sample size. Small values can be considered as conservative. Default is 0.02.
#'
#'@param blanks
#'Character vector of names of blanks. If specified, all substances found in any of the blanks
#'will be removed from all samples (i.e. c("blank1", "blank2")). The names have to correspond
#'to a name given in the first line of \code{datafile}.
#'
#'@param delete_single_peak
#'logical, determines whether substances that occur in just one sample are
#' removed or not. By default single substances are retained in chromatograms.
#'
#'@param n_iter
#' integer indicating the iterations of the core alignment algorithm. Additional replications of
#' the alignment and merging steps can be helpful to clean-up chromatograms, that otherwise show
#' some remaining peak outliers mapped to the wrong mean retention time. Inspect alignment visually
#' with a Heatmap \code{\link{gc_heatmap}}.
#'
#' @param merge_rare_peaks
#' logical determining whether peaks that are redundant by means of a similar retention times
#' are merged. If \code{TRUE} peaks are merged as long as only 5 % of the samples contain two peaks.
#' Always the peak with the higher abundance (i.e. peak area or peak height) is retained.
#'
#' @inheritParams shared_peaks
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
#'@import magrittr
#'
#' @examples
#' data(gc_peak_data)
#' gc_peak_data <- gc_peak_data[1:4]
#' out <- align_chromatograms(gc_peak_data, conc_col_name = "area", rt_col_name = "RT",
#'        rt_cutoff_low = 5, rt_cutoff_high = 45, reference = "ind3",
#'          max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03,
#'          blanks = NULL, delete_single_peak = TRUE)
#'
#'
#' @export
#'

align_chromatograms <- function(data, sep = "\t",conc_col_name=NULL, rt_col_name = NULL, write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
                                max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = NULL,
                                delete_single_peak = FALSE,n_iter=1,merge_rare_peaks=FALSE,error=error) {

###
    ### Allow the usage of another type of reference (i.e a standard that is not a sample)
    ### Possibility:
    ### One reference is called "reference", used for shifting and then eliminated
    ### Consider to allow a small deviation (e.g 0.01) for shared peaks
    ### Allow this case in form of an if-else statement
###

    ###################################################################################
    # Show the start of the alignment process and the corresponding time to the console
    ###################################################################################
    cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%H:%M:%S")),'\n\n'))


    ###############################################
    # Check if all required arguments are specified
    ################################################
    if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
    if (is.null(conc_col_name)) stop("Column containing concentration of peaks is not specified. Define conc_col_name")
    if (is.null(reference)) stop("Reference is missing. Specify a reference to align the others to")

    # check if data is a .txt file to load
    if (is.character(data)) {
        if (stringr::str_detect(data, ".txt")) {
            ###############
            # extract names
            ###############
            ind_names <- readr::read_lines(data, n_max = 1) %>%
                stringr::str_split(pattern = sep) %>%
                unlist()
            ind_names <- ind_names[ind_names != ""]    #.[. != ""]

            ######################
            # extract column names
            ######################
            col_names <- readr::read_lines(data, n_max = 1, skip = 1) %>%
                stringr::str_split(pattern = sep) %>%
                unlist()
            col_names <- col_names[col_names != ""]    #.[. != ""]
            ########################################
            # remove leading and tailing whitespaces
            ########################################
            col_names <- stringr::str_trim(col_names)
            ind_names <- stringr::str_trim(ind_names)

            ##############
            # extract data
            ##############
            gc_data <- read.table(data, skip = 2, sep = sep, stringsAsFactors = F)

            #####################
            # remove pure NA rows
            #####################
            gc_data <- gc_data[!(rowSums(is.na(gc_data)) == nrow(gc_data)), ]

            ######################
            # transform to numeric
            ######################
            gc_data <-  as.data.frame(apply(gc_data, 2, as.numeric))

            #########################################
            # Check input for completeness and format
            #########################################
            # check 1
            if (!((ncol(gc_data) / length(col_names)) %% 1) == 0) stop("Number of data columns is not a multiple of the column names provided")
            # check 2
            if (!((ncol(gc_data) / length(col_names))  == length(ind_names))) stop("Number of sample names provided does not fit to the number of columns in the data")
            # check 3
            if(any(duplicated(ind_names)))warning("Avoid duplications in sample names")

            # matrix to list
            gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)
        }

    } else if (is.list(data)) {
        # check if data is list of data.frames
        # check for data.frames
        if (!(all(unlist(lapply(data, is.data.frame))))) stop("Data object has to be a list, whereby each element is a data.frame with the GC peak data for an individual")
        # check whether all data.frames have a name
        if ((is.null(names(data)))) stop("Data object has to be a list, whereby each element has to be named with the ID of the respective individual")
        # check whether all data.frames contain the same column names
        col_names <- unlist(lapply(data, function(x) out <- names(x)))
        if(any(table(col_names) != length(data))) stop("Each data.frame in the list has to have the same variable names (i.e. 'RT' 'area')")

        ind_names <- names(data)
        col_names <- names(data[[1]])
        gc_peak_list <- data

    }

cat(paste0('GC-data for ',as.character(length(ind_names)),' samples loaded\n ',
               'Range of Relative Variation: ',as.character(round(align_var(gc_peak_list,rt_col_name)$range[1],2)),
               '\u002d',as.character(round(align_var(gc_peak_list,rt_col_name)$range[2],2))," ... Average Relative Variation: ",
               as.character(round(align_var(gc_peak_list,rt_col_name)$average,2)),
               '\n','\n'))

    ########################
    # Start of processing
    ########################

    # 1.) cut retention times
    gc_peak_list <- lapply(gc_peak_list, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_col_name)

    #####################
    # save for later use
    ####################

    gc_peak_list_raw <- lapply(gc_peak_list, matrix_append, gc_peak_list)

    if(!is.null(rt_cutoff_low) & is.null(rt_cutoff_high)){
    cat(paste0('Retention Time Cut-Off applied:\n', 'Everything below ',as.character(rt_cutoff_low),' minutes deleted','\n'))
    }

    if(!is.null(rt_cutoff_high) & is.null(rt_cutoff_low)){
        cat(paste0('Retention Time Cut-Off applied:\n', 'Everything above ',as.character(rt_cutoff_high),' minutes deleted','\n'))
    }

    if(!is.null(rt_cutoff_high) & !is.null(rt_cutoff_low)){
        cat(paste0('Retention Time Cut-Off applied:\n', 'Everything below ',as.character(rt_cutoff_low),' and above ',as.character(rt_cutoff_high) ,' minutes deleted','\n'))
    }


    # 2.) Linear Transformation of Retentiontimes

    cat('\nLinear Transformation ... ')

    ## thinking about reference: default is chromatogram with most peaks - optional: manual  !Currently it is manual
    if(reference=="reference"){ # New option
    gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift=max_linear_shift, step_size=0.01,
                                            error=error, reference = reference, rt_col_name = rt_col_name)
    gc_peak_list_linear <- gc_peak_list_linear[-which(names(gc_peak_list_linear)==reference)]; # without elements reference
    }else{
        gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift=max_linear_shift, step_size=0.01,
                                                     error=error, reference = reference, rt_col_name = rt_col_name)
        }

     cat(paste('Done ... ','\n','Range of Relative Variation: ',as.character(round(align_var(gc_peak_list_linear,rt_col_name)$range[1],2)),
              '\u002d',as.character(round(align_var(gc_peak_list_linear,rt_col_name)$range[2],2))," ... Average Relative Variation: ",
              as.character(round(align_var(gc_peak_list_linear,rt_col_name)$average,2))),'\n','\n')
gc_peak_list_linear <- lapply(gc_peak_list_linear, matrix_append, gc_peak_list_linear)

    #############
    # align peaks
    #############

    cat(c('Start Alignment of Peaks ... ','This might take a while!\n','\n'))

    Fun_Fact()
    gc_peak_list_aligned <- gc_peak_list_linear
    no_peaks <- matrix(NA,nrow = n_iter,ncol = 1)
    merged_peaks <- matrix(NA, nrow = n_iter,ncol = 1)
    for (R in 1:n_iter){ # Allows to iteratively execute the algorithm

        gc_peak_list_aligned <- align_individual_peaks(gc_peak_list_aligned, max_diff_peak2mean = max_diff_peak2mean, n_iter = n_iter, rt_col_name = rt_col_name,R=R)

    # Check whether zero rows are present
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    # delete empty rows (if existing)
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) {
        keep_rows <- which(!is.na(average_rts))
        out <- x[keep_rows, ]
    })

    # calculate average rts again for merging
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    no_peaks[R] <- nrow(gc_peak_list_aligned[[1]]) # No before merging
    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name)
    merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])
    no_peaks[R] <- nrow(gc_peak_list_aligned[[1]]) # No after merging

    if(R==n_iter & merge_rare_peaks==TRUE){
        average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

        gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name,criterion="proportional",conc_col_name = conc_col_name)
        merged_peaks[R] <- no_peaks[R] - nrow(gc_peak_list_aligned[[1]])
        rare_peak_pairs <- 0 # cause they were eliminated
#     }else if(R==n_iter){
#         average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name)
#         similar <- similar_peaks(average_rts, min_diff_peak2peak)    # remaining similarities
#         rare_peak_pairs <- length(similar) # WARNING: These are not definitive redundant, needs some thinking
    }

    cat('\n','Range of Relative Variation: ',as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$range[1],2)),
        '\u002d',as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$range[2],2))," ... Average Relative Variation: ",
        as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$average,2)),
        '\n\n')
    cat('Merged Redundant Peaks')

    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    # delete empty rows again

    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)

    } # End of iterative alignment and merging

    cat(paste('\n','Peak Alignment Done'),'\n\n')

    ###############################################
    # Sort chromatograms back to the initial order
    ###############################################
    if(reference=="reference"){
    gc_peak_list_raw <- gc_peak_list_raw[-which(names(gc_peak_list_raw)==reference)]; # remove the reference
}
    gc_peak_list_aligned <- gc_peak_list_aligned[match(names(gc_peak_list_raw),names(gc_peak_list_aligned))]



    # delete blanks
    if (!is.null(blanks)) {
        # delete one blank
        delete_blank <- function(blank, gc_peak_list_aligned) {
            del_substances <- which(gc_peak_list_aligned[[blank]][[rt_col_name]] > 0)
            chroma_out <- lapply(gc_peak_list_aligned, function(x) x[-del_substances,])
            chroma_out
        }
        # delete all blanks
        for (i in blanks) {
            gc_peak_list_aligned <- delete_blank(i, gc_peak_list_aligned)
        }
        cat('Blank Peaks deleted & Blanks removed\n\n')
    }


    # delete single substances
    # create matrix with all retention times
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))

    singular_peaks <- nrow(rt_mat) # To determine number of deleted peaks
    if (delete_single_peak) {
        # find single retention times in rows
        single_subs_ind <- which(rowSums(rt_mat > 0) == 1)

        if (length(single_subs_ind) > 0) {
            # delete substances occuring in just one individual
            gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) x[-single_subs_ind, ])
        }
        cat(paste('Single Peaks deleted:',as.character(length(single_subs_ind)),'have been removed\n\n'))

    }

    # calculate final retention times
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))
    singular_peaks <- singular_peaks - nrow(rt_mat) # how many were deleted, if any
     #######################
     # Outputs for Heatmaps
     #######################

    rt_raw <- rt_extract(gc_peak_list = gc_peak_list_raw,rt_col_name =rt_col_name)
    rt_linear <- rt_extract(gc_peak_list = gc_peak_list_linear,rt_col_name =rt_col_name)
    rt_aligned <- rt_extract(gc_peak_list = gc_peak_list_aligned,rt_col_name =rt_col_name)

    # ============================
    # mean per row without 0
    row_mean <- function(x) {
        if (all(x==0)) 0 else mean(x[x!=0])
    }
    mean_per_row <- apply(rt_mat,1, row_mean)

    # create output matrices for all variables
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


    if (!is.null(write_output)){

        ####################################################################
        # Extract the name of the datafile for proper naming of output files
        ####################################################################
        if(is.character(data)){
        prefix <- strsplit(data,split = "/")
        prefix <- as.character(prefix[[1]][length(prefix[[1]])])
        prefix <- as.character(strsplit(prefix,split = ".txt"))
        } else {
            prefix <- "Aligned"
        }
        write_files <- function(x) {
            write.table(output[[x]], file = paste0(prefix,"_", x, ".txt"), sep = "\t", row.names = FALSE)
        }
        lapply(write_output, write_files)
    }

    ### Write paramters choosen for the algorithm to a file for documentation

    align_summary <- list(No_of_samples = length(gc_peak_list_raw),
                          No_Peaks_aligned=ncol(rt_aligned)-1,
                          variance_before=align_var(gc_peak_list_raw,rt_col_name),
                          variance_aligned=align_var(gc_peak_list_aligned,rt_col_name),
                          singular_peaks=singular_peaks)

    output_algorithm <- list(aligned=output,summary=align_summary,
                             heatmap_input=list(initial_rt=rt_raw,
                                                linear_shifted_rt=rt_linear,aligned_rt=rt_aligned),
                             call=match.call())



    class(output_algorithm) <- "GCalign" # name of list


    cat(paste('Alignment was Successful!\n','Time:'),strftime(Sys.time(),format = "%H:%M:%S"),'\n')
    return(output_algorithm)
}
