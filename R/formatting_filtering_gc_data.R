add_linshifts <- function(mat = NULL, object = NULL) {
    samples <- colnames(mat)
    df <- object[["Logfile"]][["LinearShift"]]
    for (x in samples) {
mat[,which(colnames(mat) == x)][mat[,which(colnames(mat) == x)] > 0] <-
            mat[,which(colnames(mat) == x)][mat[,which(colnames(mat) == x)] > 0] +
    df[["shift"]][which(df[["sample"]] == x)]
    }
    return(mat)
}#add_linshift


add_linshifts2 <- function(dx = NULL, rt_col_name = NULL, Logbook = NULL) {
    df <- Logbook[["LinearShift"]]
    samples <- names(dx)
    for (x in samples) {
        dx[[x]][[rt_col_name]][dx[[x]][[rt_col_name]] > 0] <- dx[[x]][[rt_col_name]][dx[[x]][[rt_col_name]] > 0] + df[["shift"]][which(df[["sample"]] == x)]
    }
    return(dx)
}#add_linshift2

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
}#align_var

check_redundancy <- function(gc_peak_df, similar_peaks, rt_col_name){
    # If only one of two neighbouring rows contain a substance
    # they are redundant, coded with 1
    row1 <- gc_peak_df[similar_peaks - 1, rt_col_name]
    row2 <- gc_peak_df[similar_peaks, rt_col_name]
    redundant <- 0
    if (row1 == 0 | row2 == 0) {
        redundant <- 1
    }
    return(redundant)
}#check_redundancy

conv_gc_mat_to_list <- function(gc_data, ind_names, var_names) {
    extract <- seq(from = 1, to = ncol(gc_data), by = length(var_names))
    chromatograms <- lapply(extract, function(x) gc_data[, x:(x + length(var_names) - 1)])
    names(chromatograms) <- ind_names
    # Inserted to ensure each sample is (converted) to a data frame
    chromatograms <- lapply(X = chromatograms, as.data.frame)
    chromatograms <- lapply(chromatograms, rename_cols, var_names)
    return(chromatograms)
}#conv_gc_mat_to_list

correct_colnames <- function(gc_peak_df,col_names) {
    colnames(gc_peak_df) <- col_names
    return(gc_peak_df)
}#correct_colnames

delete_blank <- function(blanks, gc_peak_list_aligned, rt_col_name) {

    # indices of peaks
    delete <- sort(unique(unlist(lapply(blanks, function(fx) {
        which(gc_peak_list_aligned[[fx]][[rt_col_name]] > 0)
    }))))

    chroma_out <- gc_peak_list_aligned

    # remove blanks
    chroma_out[blanks] <- NULL

    # remove peaks from samples
    if (length(delete) > 0) chroma_out <- lapply(chroma_out, function(x) x[-delete,])

    return(chroma_out)
}

delete_empty_rows <- function(gc_peak_df, average_rts){
    gc_peak_df <- gc_peak_df[!is.na(average_rts), ]
    gc_peak_df
}#delete_empty_rows

delete_space_colnames <- function(gc_data) {
    names(gc_data) <- stringr::str_replace_all(names(gc_data), " ", "")
    gc_data
}#delete_space_colnames

dummy_col <- function(x) {
    x2 <- rep(NA, nrow(x))
    x2[1:length(x[!is.na(x)])] <- 1:length(x[!is.na(x)])
    x <- data.frame(x, GCalignR_Dummy = x2)
    return(x)
}#dummy_col

dummy_remove <- function(x) {
    return(x[,-2])
}#dummy_remove

function_call <- function(call,FUN="align_chromatograms"){
    form <- formals(FUN) # all arguemnts and there defaults, blank means no default.
    for (n in names(form)) { # for every args of the function
        if (!(n %in% names(call))) call <- append(call,form[n])  ## add missing args to list call
    }
    return(call)
}#function_call

is_redundant <- function(redundant, criterion="strict"){
    # Indicates by a binary output variable (1/0) if rows should be merged
    # Methods: Strict: A single sample with two peaks prevents merging
    #           Proportional: Merging is acceptabel if only 5 % of samples show two peaks
    ToMerge <- 0
    if (criterion == "strict") {
        if (sum(redundant)/length(redundant) == 1) {
            ToMerge <- 1
        }
    } else if (criterion == "proportional") {
        if (sum(redundant)/length(redundant) >= 0.95) {
            ToMerge <- 1
        }
    }
    ToMerge
}#is_redundant


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
}#matrix_append

mean_retention_times <- function(gc_peak_list, rt_col_name) {
    n_substance <- nrow(gc_peak_list[[1]])
    out <- unlist(lapply(1:n_substance,
                         function(x) mean_retention_time_row(gc_peak_list, 1:length(gc_peak_list), x, rt_col_name)))
    return(out)
}#mean_retention_times

mean_retention_time_row <- function(gc_peak_list, samples, retention_row, rt_col_name){
    xy <- function(x, retention_row, rt_col_name) {
        out <- x[retention_row, rt_col_name]
        return(out)
    }
    rts <- unlist(lapply(gc_peak_list[samples], xy,retention_row, rt_col_name))
    mean_rt <- mean(rts[!(rts == 0)], na.rm = TRUE)
    return(mean_rt
    )
}#mean_retention_time_row

merge_redundant_peaks <- function(gc_peak_list,min_diff_peak2peak=0.05, rt_col_name, conc_col_name = NULL, criterion="strict"){
    merging <- 'start'
    while (merging != 'stop') {

        # calculate mean retention times
        average_rts <- mean_retention_times(gc_peak_list, rt_col_name)
        # update similarity assessment
        similar <- similar_peaks(average_rts, min_diff_peak2peak)
        counter <- 1

        while (counter != 'stop') {
            total <- ifelse(length(similar) > 0, length(similar), 1)
            # create progress bar

            if (interactive()) {
            pb <- utils::txtProgressBar(min = 0, max = total, style = 3, char = "+", width = 80)
            utils::setTxtProgressBar(pb, ifelse(is.numeric(counter),counter, total))
            }
            # stop when there are no redundancies
            if (length(similar) == 0) {
                merging <- "stop"
                break
            }
            redundant <- sapply(lapply(gc_peak_list, check_redundancy,similar_peaks = similar[counter], rt_col_name = rt_col_name), as.vector)

            to_merge <- is_redundant(redundant = redundant, criterion = criterion)
            # prove if rows are mergeable
            if (to_merge == 1) {
                gc_peak_list <- lapply(gc_peak_list, merge_rows, to_merge = similar[counter], criterion, rt_col_name,conc_col_name)
                counter <- 'stop'
            } else if  (to_merge == 0) {
                counter <- counter + 1
                if (counter > length(similar)) {
                    merging <- 'stop'
                    counter <- 'stop'
                }
            }
        }
    }
    if (exists("pb")) close(pb)
    return(gc_peak_list)
}#merge_redundant_peaks

merge_rows <- function(gc_peak_df, to_merge, criterion="strict", rt_col_name,conc_col_name){
    # Check always the row containing just zeros, in case of zeros in both, just delete one of them
    # To Merge == Last row of a similar pair
    Row1 <- to_merge - 1
    Row2 <- to_merge
    R1 <- gc_peak_df[Row1, rt_col_name]
    R2 <- gc_peak_df[Row2, rt_col_name]
    if (criterion == "strict") {
        if (R1 == 0) {
            #  Delete Row1, if no peak exists
            gc_peak_df <- gc_peak_df[-Row1,]
        } else if (R2 == 0) {
            # Delete Row2
            gc_peak_df <- gc_peak_df[-Row2,]
        }
    }

    if (criterion == "proportional") {
        # keep the peak with higher area
        if (gc_peak_df[Row1,conc_col_name] >= gc_peak_df[Row2,conc_col_name]) {
            gc_peak_df <- gc_peak_df[-Row2,]
        } else if (gc_peak_df[Row1,conc_col_name] < gc_peak_df[Row2,conc_col_name]) {
            gc_peak_df <- gc_peak_df[-Row1,]
        }
    }
    return(gc_peak_df)
}#merge_rows

peak_counter <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for (i in 1:length(gc_peak_list)) {
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt)) & rt != 0])
    return(length(rt))
}#peak_counter

peak_lister <- function(gc_peak_list,rt_col_name){
    rt <- numeric(0)
    for (i in 1:length(gc_peak_list)) {
        rt <- c(rt,gc_peak_list[[i]][[rt_col_name]])
    }
    rt <- unique(rt[!(is.na(rt)) & rt != 0])
    return(rt)
}#peak_lister

peaks2chroma <- function(data = NULL, sample = NULL, x = seq(from = 0, to = 30, length = 10000)) {
    scale <- max(data[["y"]])
    data <- data[data[["sample"]] == sample,]
    y <- rep(0, length(x))
    for (i in 1:nrow(data)) {
        y <- y + dnorm(x, mean = data[["x"]][i], sd = 1.1 - (data[["y"]][i]/scale))
    }
    return(y)
}#peaks2chroma

remove_linshifts <- function(dx = NULL, rt_col_name = NULL, Logbook = NULL) {
    df <- Logbook[["LinearShift"]]
    samples <- names(dx[[rt_col_name]])[2:length(names(dx[[rt_col_name]]))]
        for (x in samples) {
    dx[[rt_col_name]][[x]][dx[[rt_col_name]][[x]] > 0] <- dx[[rt_col_name]][[x]][dx[[rt_col_name]][[x]] > 0] - df[["shift"]][which(df[["sample"]] == x)]
        }
    return(dx)
}#remove_linshifts

remove_linshifts2 <- function(dx = NULL, rt_col_name = NULL, Logbook = NULL) {
    df <- Logbook[["LinearShift"]]
    samples <- names(dx)
    for (x in samples) {
        dx[[x]][[rt_col_name]][dx[[x]][[rt_col_name]] > 0] <- dx[[x]][[rt_col_name]][dx[[x]][[rt_col_name]] > 0] - df[["shift"]][which(df[["sample"]] == x)]
    }
    return(dx)
}#remove_linshifts2

rename_cols <- function(data, var_names) {
    names(data) <- var_names
    data
}#rename_cols

remove_gaps <- function(gc_peak_list, rt_col_name) {
    gc_peak_list <- lapply(gc_peak_list, FUN = function(x) {
        if (any(is.na(x[[rt_col_name]]))) {
            p <- as.vector(which(is.na(x[[rt_col_name]])))
            x <- x[-p,]
        }
        if (any(x[[rt_col_name]] == 0)) {
            p <- as.vector(which(x[[rt_col_name]] == 0))
            x <- x[-p]
        }
        return(x)
    })
    return(gc_peak_list)
}#remove_gaps

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
}#rt_cutoff

rt_extract <- function(gc_peak_list, rt_col_name) {
    # blanks and del_single_sub are removed, since their removal
    # is of importance only for the last step, where it is applied
    # outside this function call
    gc_peak_list <- lapply(gc_peak_list, matrix_append, gc_peak_list)
    id <- names(gc_peak_list)
    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))
    rt_mat <- do.call(cbind, lapply(gc_peak_list, function(x) x[[rt_col_name]]))
    rt_mat <- as.data.frame(t(rt_mat))
    rt_mat2 <- rt_mat
    rt_mat2[rt_mat2 == 0] <- NA
    colnames(rt_mat) <-
    as.character(colMeans(rt_mat2,na.rm = T)) # No rounding, are not plotted as labels anyway
    rt_mat <- cbind(id,rt_mat)
}#rt_extract

rt_min_max <- function(df, rt_col_name) data.frame(min = min(df[[rt_col_name]][df[[rt_col_name]] > 0], na.rm = T), max = max(df[[rt_col_name]], na.rm = T))

conc_max <- function(df, conc_col_name) data.frame(min = min(df[[conc_col_name]][df[[conc_col_name]] > 0]), max = max(df[[conc_col_name]]))

p2c <- function(df = NULL, x = NULL, conc_max = NULL, rt_col_name = NULL, conc_col_name = NULL, width = NULL) {
    temp <- rep(0, length(x))
    if (is.null(conc_max)) for (i in 1:nrow(df)) temp <- temp + dnorm(x, mean = df[[rt_col_name]][i], sd = width)
    if (!is.null(conc_max)) for (i in 1:nrow(df)) temp <- temp + dnorm(x, mean = df[[rt_col_name]][i], 1.1 - df[[conc_col_name]][i]/conc_max)
    return(temp)
}#rt_min_max

shift_rows = function(chromatograms, current_sample_index, retention_row) {
    n_col <- ncol(chromatograms[[1]])
    zeros <- as.data.frame(matrix(0,nrow = 1,ncol = n_col))
    colnames(zeros) <- names(chromatograms[[1]])
    chroma_temp <-  chromatograms[[current_sample_index]]

    if (retention_row != 1) {
        chroma_temp <- rbind(chroma_temp[1:(retention_row - 1),], zeros,
                             chroma_temp[retention_row:nrow(chroma_temp), ])
    } else {
        chroma_temp <- rbind(zeros, chroma_temp)
    }

    chromatograms[[current_sample_index]] <- chroma_temp
    return(chromatograms)
}#shift_rows

similar_peaks <- function(average_rts, min_diff_peak2peak = 0.05){
    difference <- rep(NA, (length(average_rts) - 1))
    for (i in 2:length(average_rts)) {
        difference[i] <- average_rts[i] - average_rts[i - 1]
    }
    # Find rows that show similar mean retention times
    similar <- which(difference <= min_diff_peak2peak)
    return(similar)
}#similar_peaks

time_cut <- function(df, rt_col_name, rt_limits) {
    min <- min(rt_limits)
    max <- max(rt_limits)
    df <- subset(df, (df[[rt_col_name]] >= min) & (df[[rt_col_name]] <= max))
    return(df)
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
}#write_files
