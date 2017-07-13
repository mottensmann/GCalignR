context("align_peaks")

library(testthat)
data("peak_data")
data <- peak_data[1:3]

set.seed(1)
x <- GCalignR:::align_peaks(gc_peak_list = data, rt_col_name = "time")

test_that("output is correct", {
})
