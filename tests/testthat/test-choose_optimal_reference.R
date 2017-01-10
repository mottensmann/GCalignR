context("choose_optimal_reference")

library(testthat)
data("peak_data")

out <- choose_optimal_reference(gc_peak_list = peak_data[1:3],rt_col_name = "time")
out2 <- choose_optimal_reference(gc_peak_list = peak_data[1:4],rt_col_name = "time")

test_that("output is correct", {
    expect_equal(out[["sample"]], "C3") # picked reference
    expect_equal(out[["score"]], 56) # Median number of shared peaks
    expect_equal(out2[["sample"]], "C3") # picked reference
    expect_equal(out2[["score"]], 40.5) # Median number of shared peaks
})







