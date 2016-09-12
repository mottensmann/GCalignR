context("check_input")

library(testthat)
data("peak_data")
x <- check_input(peak_data,reference = NULL,list_peaks = TRUE,main = "TEST", xlab = "Samples",  ylab = "Peaks")

test_that("output is correct", {
    expect_equal(as.character(x[["sample"]][[78]]), "P41")
    expect_equal(x[["peaks"]][28], 117)
})




