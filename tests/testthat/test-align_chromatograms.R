context("align_chromatograms")

library(testthat)
data("gc_peak_data")
data <- gc_peak_data[1:3]

set.seed(1)
out <- align_chromatograms(data, sep = "\t", conc_col_name= "area", rt_col_name = "RT", write_output = NULL,
        rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "ind1",
        max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = NULL,
        delete_single_peak = FALSE, n_iter=1)

test_that("output is correct", {
    expect_equal(out$summary$No_Peaks_aligned, 162)
    expect_equal(round(out$summary$variance_aligned$average, 6), 0.000909)
})







