context("linear_transformation")

library(testthat)
data("peak_data")
data <- peak_data[1:5]
data <- lapply(data, function(x) x[1:10,])

set.seed(1)
out <- align_chromatograms(data, sep = "\t", rt_col_name = "time", write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "M2",max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = "C2",delete_single_peak = FALSE)

set.seed(1)
out2 <- align_chromatograms(data, sep = "\t", rt_col_name = "time", write_output = NULL, rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "M2",max_linear_shift = 0, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = "C2",delete_single_peak = FALSE)

test_that("output is correct", {
    expect_equal(sum(out[["heatmap_input"]][["linear_transformed_rts"]][,2:10]), 214.86) #214.95
    expect_equal(sum(out[["heatmap_input"]][["input_rts"]][,2:10]), 214.41) #291
    expect_equal(sum(out2[["heatmap_input"]][["linear_transformed_rts"]][,2:10]), 214.41)
    expect_equal(sum(out2[["heatmap_input"]][["input_rts"]][,2:10]), 214.41) #291
})
