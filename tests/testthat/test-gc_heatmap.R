context("gc_heatmap")

library(testthat)
data("seal_peaks_aligned")
x <- seal_peaks_aligned
out <- gc_heatmap(x)

test_that("output is correct", {
    expect_equal(length(out[["data"]][["substance"]]), 18984)

})




