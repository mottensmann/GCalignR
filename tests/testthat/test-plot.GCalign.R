context("plot.GCalign")

library(testthat)
data("aligned_peak_data")
x <- aligned_peak_data
out1 <- plot(x)
out2 <- plot(x,which_plot = 'linear_shifts')
out3 <- plot(x,which_plot = 'rt_range')
out4 <- plot(x,which_plot = 'linear_shifts')
out5 <- plot(x,which_plot = 'linear_shifts')


test_that("output is correct", {
    expect_equal(out1, 1)
    expect_equal(sd(out2[["shifts"]]), 0.005371991)
    expect_equal(out3[3,1], 0.04)

})




