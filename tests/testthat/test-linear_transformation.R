context("linear_transformation")

library(testthat)
data("peak_data")
data <- peak_data[1:3]

set.seed(1)
out <- align_chromatograms(data, sep = "\t", rt_col_name = "time", write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "M2",max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = "C2",delete_single_peak = FALSE)


test_that("output is correct", {
    expect_equal(out[["aligned"]][["time"]][["C3"]][6], 5.52)
    expect_equal(out[["aligned"]][["time"]][["C3"]][53], 15.17)
})
