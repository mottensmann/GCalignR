#' Clean-Up and formatting of Peak lists
#'
#' @description
#' Set of functions used internally to to format the data for efficient processing.
#' Steps include the removal of redundant rows (i.e. no sample carries a substance)
#'
#' @keywords internal
#'
#'
align_var <- function(gc_peak_list,rt_col_name){
    # Calculates the range of retention times for each peak, estimate is the width
    # computated as the distance between min and max values
    width_of_peaks <- function(gc_peak_list, sample_indices, peak, rt_col_name){
        rt <- unlist(lapply(gc_peak_list[sample_indices], function(x) x[peak, rt_col_name])) # all rts per peak
        if (any(!is.na(rt))) { # Check if all are empty
            min2max <- range(rt[!(rt == 0)], na.rm = TRUE) # range of rts, width per peak
            width <- abs(diff(min2max)) # If zero, no deviation, or just one sample!
        } else {
            width <- NA
        }
        return(width)
    }
    peak <- nrow(gc_peak_list[[1]]) # number of peaks
    out <- unlist(lapply(1:peak,
                         function(x) width_of_peaks(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    output <- list(range = round(range(out,na.rm = T),2),Average = round(mean(out,na.rm = T),2), Std = round(stats::sd(out,na.rm = T),2))
    output <- unlist(output)
    names(output)[1:2] <- c("Min","Max")
    return(output)
}

conv_gc_mat_to_list <- function(gc_data, ind_names, var_names) {
    extract <- seq(from = 1, to = ncol(gc_data), by = length(var_names))
    chromatograms <- lapply(extract, function(x) gc_data[, x:(x + length(var_names) - 1)])
    names(chromatograms) <- ind_names
    # Inserted to ensure each sample is (converted) to a data frame
    chromatograms <- lapply(X = chromatograms, as.data.frame)
    chromatograms <- lapply(chromatograms, rename_cols, var_names)
    return(chromatograms)
}

correct_colnames <- function(gc_peak_df,col_names) {
    colnames(gc_peak_df) <- col_names
    return(gc_peak_df)
}#end

delete_empty_rows <- function(gc_peak_df, average_rts){
    gc_peak_df <- gc_peak_df[!is.na(average_rts), ]
    gc_peak_df
}

delete_space_colnames <- function(gc_data) {
    names(gc_data) <- stringr::str_replace_all(names(gc_data), " ", "")
    gc_data
}

dummy_col <- function(x) {
    x2 <- rep(NA, nrow(x))
    x2[1:length(x[!is.na(x)])] <- 1:length(x[!is.na(x)])
    x <- data.frame(x, GCalignR_Dummy = x2)
    return(x)
}

dummy_remove <- function(x) {
    return(x[,-2])
}
function_call <- function(call,FUN="align_chromatograms"){
    form <- formals(FUN) # all arguemnts and there defaults, blank means no default.
    for (n in names(form)) { # for every args of the function
        if (!(n %in% names(call))) call <- append(call,form[n])  ## add missing args to list call
    }
    # type <- as.vector(which(lapply(call, function(x) out <- class(x))!="NULL"))
    # call[type] <- lapply(call[type], function(x) x <- as.character(x)) # if not NULL, convert to char
    # call[-type] <- lapply(call[-type], function(x) x <- "NULL") # if NULL --> "NULL"
    # call <- do.call(rbind,call) # creates data frame > 1 column if more than one blank!
    # call <- t(as.data.frame(call))
    # row.names(call) <- NULL
    return(call)
}

matrix_append <- function(gc_peak_df, gc_peak_list,val = c("Zero","NA")) {
    val <- match.arg(val)
    val <- ifelse(val == "Zero",0,NA)
    # Add zeros or NAs to matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(gc_peak_list,function(x) nrow(x)))
    ToAppend <- MaxLength - nrow(gc_peak_df)
    Cols <- ncol(gc_peak_df)
    Zeros <- matrix(val,nrow = ToAppend,ncol = Cols)
    colnames(Zeros) <- names(gc_peak_df)
    gc_peak_df <- rbind(gc_peak_df[,],Zeros)
    return(gc_peak_df)
}

peak_counter <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for (i in 1:length(gc_peak_list)) {
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt)) & rt != 0])
    return(length(rt))
}

peak_lister <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for (i in 1:length(gc_peak_list)) {
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt)) & rt != 0])
    return(rt)
}

remove_linshifts <- function(dx = NULL, rt_col_name = NULL, Logbook = NULL) {
    df <- Logbook[["LinearShift"]]
    samples <- names(dx[[rt_col_name]])[2:length(names(dx[["time"]]))]
        for (x in samples) {
    dx[[rt_col_name]][[x]][dx[[rt_col_name]][[x]] > 0] <- dx[[rt_col_name]][[x]][dx[[rt_col_name]][[x]] > 0] - df[["shift"]][which(df[["sample"]] == x)]
        }
    return(dx)
}

rename_cols = function(data, var_names) {
    names(data) <- var_names
    data
}

rt_cutoff <- function(gc_peak_df, rt_col_name, low = NULL, high = NULL) {
    # make sure gc_peak_df is a data frame
    var_names <- names(gc_peak_df)
    gc_peak_df <- as.data.frame(gc_peak_df)
    gc_peak_df <- rename_cols(data = gc_peak_df, var_names = var_names)

    # last row in the dataset, i.e. maximum number of peaks per sample
    highrow <- nrow(gc_peak_df)
    lowrow <- 1
    if (!is.null(low)) {
        lowrow <- min(which(gc_peak_df[[rt_col_name]] > low))
    }
    if (!is.null(high)) {
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
    rt_mat2[rt_mat2 == 0] <- NA
    colnames(rt_mat) <-
        as.character(colMeans(rt_mat2,na.rm = T)) # No rounding, are not plotted as labels anyway
    # colnames(rt_mat) <- as.character(1:ncol(rt_mat)) # does not work cause gc_heatmap needs numbers!

    rt_mat <- cbind(id,rt_mat)
}

write_files <- function(var = NULL, data = NULL, name = NULL) {
    filename <- paste0(name,"_",var, ".txt")
    c <- 1
    while (file.exists(filename)) {
        filename <- paste0(name,"_",var,"_",as.character(c),".txt")
        c <- c + 1
    }
    utils::write.table(data[[var]], # change to [[]]
                       file = filename,
                       sep = "\t",
                       row.names = FALSE)
    return(filename)
}

