


align_chromatograms <- function(data, ind_names, var_names = c("RT", "startRT", "endRT", "Area", "xArea", "Height", "xHeight"),
                                rt_cutoff_low = 8,) {
    
    # 1.) cut retention times below 8
    chromatograms <- lapply(chromatograms, rt_cutoff, low = 11, high = 20, rt_col_name = "RT")
    

    # 2.) Linear Transformation of Retentiontimes
    chroma_aligned <- linear_transformation(chromatograms, shift=0.05, step_size=0.01, error=0, reference = "w3", rt_col_name = "RT")
    
    # Make List equal in length
    chromatograms <- lapply(chroma_aligned, matrix_append, chroma_aligned)
    

    Length <- (max(unlist(lapply(chromatograms, function(x) out <- nrow(x))))) # To obtain Rows after run of the algorithm
    Variation <- mean(var_per_row(chromatograms),na.rm = T)
    
    chromatograms_aligned <- align_individual_peaks(chromatograms, error_span = 0.02, n_iter = 1)
    
}