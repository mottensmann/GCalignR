#' Aligning chromatograms based on retention times
#'
#'@description
#'\code{align_chromatograms()} is the core function of \link{GcalignR}.
#'
#'@details
#'Alignment is archieved by running three major algorithms always considering the complete
#'set of samples submitted to the function. In brief: All chromatograms are linearly shifted to
#'maximise similarity with a common reference to account for systematic shifts in retention times
#'caused by the gas-chromatography setup. Second peaks of similar retention times are step-wise
#'moved to the same location (i.e row) to group substances. In a third step peaks (i.e.) substances
#'that show smaller differences in mean retention times than expected by the achievable resolution
#'of gas-chromatography or the chemistry of the compounds are merged. Optional processing includes the
#'removal of peaks present in blanks (i.e. contaminations) and peaks that a uniquely found within just
#'one sample.
#'
#'@param datafile path to a datafile. It has to be a .txt file. The first row needs to contain
#'  sample names, the second row column names of the corresponding chromatograms. Starting with the
#'  third row chromatograms are included, whereby single samples are concatenated horizontally. Each
#'  chromatogram needs to consist of the same number of columns, at least two are required:The
#'  retention time and a measure of concentration (e.g. peak area or height).
#'
#'@param sep the field separator character. Values on each line of the file are separated by this
#'  character. The default is tab seperated (sep = '\\t'). See sep argument in read.table for details.
#'
#'@param conc_col_name
#'character string naming a column used for quantification of peaks (e.g. peak area or peak height).
#'The designated variable needs to be numeric.
#'
#'@param rt_col_name
#'character string naming the column holding retention times.The designated variable needs to be numeric.
#'
#'@param write_output character vector of variables to write to a text file (e.g. \code{c("RT","Area")}.
#' The default is \code{NULL}.Names of the text files are concatenations of \code{datafile} and the output variables.
#'
#'@param rt_cutoff_low lower threshold under which retention times are cutted (i.e. 5)
#'
#'@param rt_cutoff_high upper threshold above which retention times are cutted (i.e. 35)
#'
#'@param reference \code{character}, a sample(name) to which all other samples are aligned to by means of a
#'  linear shift (i.e. "individual3") The name has to correspond to an individual name given
#'  in the first line of \code{datafile}.
#'
#'@param max_linear_shift defines a window to search for an optimal linear shift of samples
#'  with respect to the reference. Shifts are evaluated within - max_linear_shift to + max_linear_shift.
#'  Default is 0.05.
#'
#'@param max_diff_peak2mean defines the allowed deviation of retention times around the mean of
#'  the corresponding row. Default is 0.02.
#'
#'@param min_diff_peak2peak defines the minimum difference in retention times among distinct
#'  substances. Substances that differ less, are merged if every sample contains either one
#'  or none of the respective compounds. This parameter is a major determinant in the classification
#'  of distinct peaks. Therefore careful consideration is required to adjust this setting to
#'  your needs (i.e. the resolution of your gas-chromatography pipeline). Default is 0.02.
#'
#'@param blanks character vector of names of blanks. If specified, all substances found in any of the blanks
#'  will be removed from all samples (i.e. c("blank1", "blank2")). The names have to correspond
#'  to a name given in the first line of \code{datafile}.
#'
#'@param delete_single_peak logical, determines whether substances that occur in just one sample are
#'  removed or not. By default single substances are retained in chromatograms.
#'
#'@param n_iter integer indicating the iterations of the core alignment algorithm.
#'
#'@return
#' Returns an object of class GCalign that is a a list with the following elements:
#' \item{call}{function call}
#' \item{chroma_aligned}{a list containing data.frames with the aligned variables}.
#' \item{rt_raw}{a data.frame with the retention times before alignment}
#' \item{rt_linear}{a data.frame with the retention times after the linear transformation}
#' \item{rt_aligned}{a data.frame with the final aligned retention times}
#' \item{parameter}{a data.frame}{containing the arguments of the function call}
#'
#'@author Martin Stoffel (martin.adam.stoffel@@gmail.com) & Meinolf Ottensmann
#'  (meinolf.ottensmann@@web.de)
#'
#'@import magrittr
#'
#'
#' @export
#'

align_chromatograms <- function(datafile, sep = "\t",conc_col_name=NULL, rt_col_name = NULL, write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = NULL,
                                max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = NULL,
                                delete_single_peak = FALSE,n_iter=1) {

    ##########################################################
    # start a clock to estimate time the function takes to run
    ##########################################################
    start.time <- pracma::tic()

    ###################################################################################
    # Show the start of the alignment process and the corresponding time to the console
    ###################################################################################
    cat(paste0('Run GCalignR\n','Start: ',as.character(strftime(Sys.time(),format = "%H:%M:%S")),'\n####################','\n','\n'))

    ###############################################
    # Check if all required arguments are specified
    ################################################
if (is.null(rt_col_name)) stop("Column containing retention times is not specifed. Define rt_col_name")
if (is.null(conc_col_name)) stop("Column containing concentration of peaks is not specified. Define conc_col_name")
if (is.null(reference)) stop("Reference is missing. Specify a reference to align the others to")

    ###############
    # extract names
    ###############
    ind_names <- readr::read_lines(datafile, n_max = 1) %>%
        stringr::str_split(pattern = sep) %>%
        unlist() %>%
        .[. != ""]

    ######################
    # extract column names
    ######################
    col_names <- readr::read_lines(datafile, n_max = 1, skip = 1) %>%
        stringr::str_split(pattern = sep) %>%
        unlist() %>%
        .[. != ""]
    ########################################
    # remove leading and tailing whitespaces
    ########################################
    col_names <- stringr::str_trim(col_names)
    ind_names <- stringr::str_trim(ind_names)

    ##############
    # extract data
    ##############
    gc_data <- read.table(datafile, skip = 2, sep = sep, stringsAsFactors = F)

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
   # if(any(duplicated(ind_names)))stop("Duplicated sample names exist. Use unique names for samples")

    # matrix to list
    gc_peak_list <- conv_gc_mat_to_list(gc_data, ind_names, var_names = col_names)

    cat(paste0('GC-data for ',as.character(length(ind_names)),' samples loaded ... ',
               'range of CV: ',as.character(round(align_var(gc_peak_list,rt_col_name)$range[1],2)),
               '\u002d',as.character(round(align_var(gc_peak_list,rt_col_name)$range[2],2))," ... average CV: ",
               as.character(round(align_var(gc_peak_list,rt_col_name)$average,2)),
               '\n################################################################################','\n','\n'))

    ########################
    # Start of processing
    ########################

    # 1.) cut retention times
    gc_peak_list <- lapply(gc_peak_list, rt_cutoff, low = rt_cutoff_low, high = rt_cutoff_high, rt_col_name = rt_col_name)

    #####################
    # save for later use
    ####################

    gc_peak_list_raw <- lapply(gc_peak_list, matrix_append, gc_peak_list)

    if(!is.null(rt_cutoff_low)){
    cat(paste0('Retention time cut-off applied:\n', 'Everything below ',as.character(rt_cutoff_low),'minutes deleted','\n##########','\n','\n'))
    }

    if(!is.null(rt_cutoff_high)){
        cat(paste0('Retention time cut-off applied:\n', 'Everything above ',as.character(rt_cutoff_high),'minutes deleted'))
    }

    if(!is.null(rt_cutoff_high) & !is.null(rt_cutoff_low)){
        cat(paste0('Retention time cut-off applied:\n', 'Everything below ',as.character(rt_cutoff_low),' and above ',as.character(rt_cutoff_high) ,' minutes deleted'))
    }


    # 2.) Linear Transformation of Retentiontimes

    cat('Linear Transformation ... ')

    ## thinking about reference: default is chromatogram with most peaks - optional: manual  !Currently it is manual
    gc_peak_list_linear <- linear_transformation(gc_peak_list, max_linear_shift=max_linear_shift, step_size=0.01,
                                            error=0, reference = reference, rt_col_name = rt_col_name)
    cat(paste('Done ... ','Range of CV: ',as.character(round(align_var(gc_peak_list_linear,rt_col_name)$range[1],2)),
              '\u002d',as.character(round(align_var(gc_peak_list_linear,rt_col_name)$range[2],2))," ... average CV: ",
              as.character(round(align_var(gc_peak_list_linear,rt_col_name)$average,2)),
              '\n'),'######################################################################################','\n')
# floor(pracma::toc(echo = F)[[1]]/60),'minutes since start',
    # Make List equal in length
    gc_peak_list_linear <- lapply(gc_peak_list_linear, matrix_append, gc_peak_list_linear)

    #############
    # align peaks
    #############

    cat(c('Start alignment of peaks ... ','This might take a while!\n',
          '=====================================================','\n'))

    # Fun_Fact() #Maybe within each iteration
    gc_peak_list_aligned <- gc_peak_list_linear
    for (R in 1:n_iter){ # several iteration of algorithm

        gc_peak_list_aligned <- align_individual_peaks(gc_peak_list_aligned, max_diff_peak2mean = max_diff_peak2mean, n_iter = n_iter, rt_col_name = rt_col_name,R=R)


    # see whether zero rows are present
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    # delete empty rows (if existing)
    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) {
        keep_rows <- which(!is.na(average_rts))
        out <- x[keep_rows, ]
    })

    # calculate average rts again for merging
    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    # merging step
    # min distance here is crucial --> depends on sample size


    gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name)
    if(R==n_iter){
        average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

        gc_peak_list_aligned <- merge_redundant_peaks(gc_peak_list_aligned, min_diff_peak2peak=min_diff_peak2peak, rt_col_name = rt_col_name,criterion="proportional",conc_col_name = conc_col_name)

    }
    cat('range of CV: ',as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$range[1],2)),
        '\u002d',as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$range[2],2))," ... average CV: ",
        as.character(round(align_var(gc_peak_list_aligned,rt_col_name)$average,2)),
        '\n')
    cat('Merged redundant peaks ... ')

    cat(paste('Elapsed time:',floor(pracma::toc(echo = F)[[1]]/60),'minutes',
        '\n##########################################################################################','\n')) #floor(pracma::toc(echo = F)[[1]]/60),'minutes since start',


    average_rts <- mean_retention_times(gc_peak_list_aligned, rt_col_name = rt_col_name)

    # delete empty rows again

    gc_peak_list_aligned <- lapply(gc_peak_list_aligned, delete_empty_rows, average_rts)

    } # Endo of iterative alignment and merging

    cat(paste('\n','Peak alignment Done ... '),
                '\n##################################################','\n','\n')

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
    }


    # delete single substances
    # create matrix with all retention times
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))

    if (delete_single_peak) {
        # find single retention times in rows
        single_subs_ind <- which(rowSums(rt_mat > 0) == 1)

        if (length(single_subs_ind) > 0) {
            # delete substances occuring in just one individual
            gc_peak_list_aligned <- lapply(gc_peak_list_aligned, function(x) x[-single_subs_ind, ])
        }

    }


    # calculate final retention times
    rt_mat <- do.call(cbind, lapply(gc_peak_list_aligned, function(x) x[[rt_col_name]]))

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

    # Extract the name of the data source to write for example plover_area.txt instead
    # of just an unprecise area.txt

    prefix <- strsplit(datafile,split = "/")
    prefix <- as.character(prefix[[1]][length(prefix[[1]])])
    prefix <- as.character(strsplit(prefix,split = ".txt"))

    if (!is.null(write_output)){
        write_files <- function(x) {
            write.table(output[[x]], file = paste0(prefix,"_", x, ".txt"), sep = "\t", row.names = FALSE)
        }
        lapply(write_output, write_files)
    }

    ### Write paramters choosen for the algorithm to a file for documentation

    parameter <- data.frame(datafile=datafile,
                            reference=reference,
                            max_linear_shift=max_linear_shift,
                            max_diff_peak2mean=max_diff_peak2mean,
                            min_diff_peak2peak=min_diff_peak2peak,
                            delete_single_peak=delete_single_peak,
                            initial_maximum_substances=ncol(rt_raw),
                            number_of_substances_final=ncol(rt_aligned)-1)
    if(!is.null(rt_cutoff_high)){
        parameter['rt_cutoff_high'] <- rt_cutoff_high
    }
    if(!is.null(rt_cutoff_low)){
        parameter['rt_cutoff_low'] <- rt_cutoff_low
    }
    if(!is.null(blanks)){
        parameter['blanks'] <- paste(blanks, collapse = '; ')
    }


    #output
    output_algorithm <- list(call=match.call(),
                chroma_aligned = output, rt_raw = rt_raw, rt_linear = rt_linear,
                rt_aligned = rt_aligned,parameter=parameter)

    class(output_algorithm) <- "GCalign" # name of list

    end.time <- pracma::toc(echo = F)
    cat(paste('Alignment was successful!\n','Time:'),strftime(Sys.time(),format = "%H:%M:%S"),'\nElapsed time:',as.character(floor(end.time[[1]]/60)),'minutes','\n')
    return(output_algorithm)
}




