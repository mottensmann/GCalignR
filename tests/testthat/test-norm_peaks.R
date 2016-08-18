context("norm_peaks")

library(testthat)
data("aligned_peak_data")

x <- norm_peaks(aligned_peak_data,conc_col_name = "area")
y <- as.vector(unlist(lapply(x,sum)))
y <- any(round(y,digits = 10)!=100)

test_that("output is correct", {
    expect_equal(y, FALSE)

})







