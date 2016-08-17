context("check_input")

library(testthat)
data("peak_data")
x <- check_input(peak_data)

test_that("output is correct", {
    expect_equal(x, "TRUE")

})




