context("align_chromatograms")


data("gc_peak_data")
data <- gc_peak_data[1:2]

out <- align_chromatograms(data, sep = "\t", conc_col_name= "area", rt_col_name = "RT", write_output = NULL,
        rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "ind1",
        max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = NULL,
        delete_single_peak = FALSE, n_iter=1)

# test_that("variance is minimized", {
#     expect_equal()
# })

out <- align_chromatograms(data, sep = "\t", conc_col_name= "area", rt_col_name = "RT", write_output = NULL,
    rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "ind1",
    max_linear_shift = 0, max_diff_peak2mean = 0, min_diff_peak2peak = 0, blanks = NULL,
    delete_single_peak = FALSE, n_iter=1)

# test_that("zeros for all three main parameters steps work", {
#     expect_equal()
# })
#

out <- align_chromatograms(data, sep = "\t", conc_col_name= "area", rt_col_name = "RT", write_output = NULL,
    rt_cutoff_low = NULL, rt_cutoff_high = NULL, reference = "ind1",
    max_linear_shift = 0.05, max_diff_peak2mean = 0.02, min_diff_peak2peak = 0.03, blanks = NULL,
    delete_single_peak = FALSE, n_iter=2)

# test_that("more than one iteration works", {
#     expect_equal()
# })




