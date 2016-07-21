#' Clean-Up and formatting of Peak lists
#'
#' @description
#' Set of functions used internally to to format the data for efficient processing.
#' Steps include the removal of redundant rows (i.e. no sample carries a substance)
#'
#' @keywords internal
#'
#'
delete_empty_rows <- function(gc_peak_df, average_rts){
    gc_peak_df <- gc_peak_df[!is.na(average_rts), ]
    gc_peak_df
}

delete_space_colnames <- function(gc_data) {
    names(gc_data) <- stringr::str_replace_all(names(gc_data), " ", "")
    gc_data
}

conv_gc_mat_to_list <- function(gc_data, ind_names, var_names) {
    extract <- seq(from = 1, to = ncol(gc_data), by = length(var_names))
    chromatograms <- lapply(extract, function(x) gc_data[, x:(x+length(var_names)-1)])
    names(chromatograms) <- ind_names
    rename_cols = function(data, var_names){
        names(data) <- var_names
        data
    }

    chromatograms <- lapply(chromatograms, rename_cols, var_names)
    return(chromatograms)
}

matrix_append <- function(gc_peak_df, gc_peak_list){
    # Add zeros matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(gc_peak_list,function(x) nrow(x)))
    ToAppend <- MaxLength-nrow(gc_peak_df)
    Cols <- ncol(gc_peak_df)
    Zeros <- matrix(0,nrow=ToAppend,ncol=Cols)
    colnames(Zeros) <- names(gc_peak_df)
    gc_peak_df<- rbind(gc_peak_df[,],Zeros)
    return(gc_peak_df)
}


rt_cutoff <- function(gc_peak_df, rt_col_name, low=NULL, high=NULL){
    # RetentionCutoff removes all Retention Times below the Threshold specified by Low (default 8s).
    # In addition Retention Times above a time defined by the Value of High (Default is Null)
    # can be applied.
    highrow <- nrow(gc_peak_df)
    lowrow <- 1
    if (!is.null(low)){
        lowrow <- min(which(gc_peak_df[[rt_col_name]] > low))
    }
    if (!is.null(high)){
        highrow <- max(which(gc_peak_df[[rt_col_name]] < high))
    }

    gc_peak_df <- gc_peak_df[lowrow:highrow, ]
    out <- gc_peak_df[!is.na(gc_peak_df[, rt_col_name]), ]
    return(out)
}

rt_extract <- function(gc_peak_list,rt_col_name){

    # blanks and del_single_sub are removed, since their removal
    # is of importance only for the last step, where it is applied
    # outside this function call

    ###############################################
    # Make length equal, if differences are present
    ###############################################

    gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)
    id <- names(gc_peak_list)

    ####################################################################
    # optional, depends on arguments regarding blanks and del_single_sub
    ####################################################################

    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))


    #################################
    # calculate final retention times
    #################################

    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))
    rt_mat <- as.data.frame(t(rt_mat))
    rt_mat2 <- rt_mat
    rt_mat2[rt_mat2==0] <- NA
    colnames(rt_mat) <-
        as.character(colMeans(rt_mat2,na.rm = T)) # No rounding, are not plotted as labels anyway
    # colnames(rt_mat) <- as.character(1:ncol(rt_mat)) # does not work cause gc_heatmap needs numbers!

    rt_mat <- cbind(id,rt_mat)
}

align_var <- function(gc_peak_list,rt_col_name){
    # Calculates the range of retention times for each peak, estimate is the width
    # computated as the distance between min and max values
    width_of_peaks <- function(gc_peak_list, sample_indices, peak, rt_col_name){
        rt <- unlist(lapply(gc_peak_list[sample_indices], function(x) x[peak, rt_col_name])) # all rts per peak
        if(any(!is.na(rt))){ # Check if all are empty
            min2max <- range(rt[!(rt == 0)], na.rm = TRUE) # range of rts, width per peak
            width <- abs(diff(min2max)) # If zero, no deviation, or just one sample!
        } else{
            width <- NA
        }
        return(width)
    }
    peak <- nrow(gc_peak_list[[1]]) # number of peaks
    out <- unlist(lapply(1:peak,
                         function(x) width_of_peaks(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    output <- list(range=round(range(out,na.rm = T),2),average=round(mean(out,na.rm=T),2),std=round(sd(out,na.rm = T),2))
    return(output)
}

peak_counter <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for(i in 1:length(gc_peak_list)){
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt))&rt!=0])
    return(length(rt))
}

peak_lister <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for(i in 1:length(gc_peak_list)){
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt))&rt!=0])
    return(rt)
}

function_call <- function(call,FUN="align_chromatograms"){
    form <- formals(FUN)
    for ( n in names(form)){ # for every args of the function
        if (!(n %in% names(call))){ # Find args not called
            call <- append(call,form[n])  ## add missing args
        }
    }
    type <- as.vector(which(lapply(call, function(x) out <- class(x))!="NULL"))
    call[type] <- lapply(call[type], function(x) x <- as.character(x))
    call[-type] <- lapply(call[-type], function(x) x <- "NULL")
    call <- do.call(rbind,call)
    call <- t(as.data.frame(call))
    row.names(call) <- NULL
    return(call)
}

