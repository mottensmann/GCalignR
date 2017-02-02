context("norm_peaks")

library(testthat)
data("aligned_peak_data")

x <- norm_peaks(aligned_peak_data, conc_col_name = "area",out = "list")
y <- as.vector(unlist(lapply(x,sum)))
y <- any(round(y,digits = 10) != 100)

x2 <- norm_peaks(aligned_peak_data, conc_col_name = "area",percent = FALSE, out = "list")
y2 <- as.vector(unlist(lapply(x2, sum)))
y2 <- any(round(y2,digits = 10) != 1)


test_that("output is correct", {
    expect_equal(y, FALSE)
    expect_equal(y2, FALSE)
})







