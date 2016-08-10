context("align_chromatograms")

library(testthat)
data("seal_peaks")
data <- seal_peaks[1:3]

set.seed(1)
out <- align_chromatograms(data, sep = "\t", conc_col_name= "area", rt_col_name = "time", write_output = NULL,
        rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "M2",
        max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = "C2",
        delete_single_peak = FALSE, n_iter=1)

test_that("output is correct", {
    expect_equal(out[["Logfile"]][["Aligned"]][["Retained"]], 114)
    expect_equal(out[["Logfile"]][["Variation"]][["Aligned"]][["Max"]], 0.01)
    expect_equal(names(out), c("aligned","heatmap_input","Logfile"))
    expect_equal(any(is.na(out[["aligned"]])),FALSE)
})







