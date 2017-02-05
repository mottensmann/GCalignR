context("norm_peaks")

library(testthat)
data("aligned_peak_data")

x <- norm_peaks(aligned_peak_data, conc_col_name = "area",out = "list")
y <- as.vector(unlist(lapply(x,sum)))
## supposed to be an artefact, values deviate from 100 by approx. e-15
y <- any(round(y,digits = 10) != 100)

x2 <- norm_peaks(aligned_peak_data, conc_col_name = "area", rt_col_name = "time", out = "data.frame")
## supposed to be an artefact, values deviate from 100 by approx. e-15
y2 <- as.vector(rowSums(x2))
y2 <- any(round(y2) != 100)

test_that("output is correct", {
    expect_equal(y, FALSE)
    expect_equal(y2, FALSE)
})







