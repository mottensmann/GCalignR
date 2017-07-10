context("draw_chromatogram")

library(testthat)
data("aligned_peak_data")

x <- draw_chromatogram(data = aligned_peak_data, rt_col_name = "time", samples = "M2", show_num = FALSE, plot = FALSE)
x2 <- draw_chromatogram(data = peak_data, rt_col_name = "time", samples = c("C2","C3"), width = 0.05, plot = F)

test_that("output is correct", {
    expect_equal(length(x), 2) # list of two elements
    expect_equal(dim(x[["df"]]), c(8622, 3))
    expect_equal(sum(x2[["df"]][["n"]]), 125) # number of peaks
    expect_equal(round(sum(x2[["df"]][["y2"]]),3), 1283.398)
})







