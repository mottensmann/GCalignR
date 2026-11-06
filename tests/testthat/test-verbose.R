context("alignment progress output")

test_that("align_chromatograms can run without progress output", {
    data("peak_data", package = "GCalignR")
    peak_subset <- lapply(peak_data[1:3], function(x) x[1:50, ])

    verbose_result <- NULL
    expect_output(
        verbose_result <- align_chromatograms(
            peak_subset,
            rt_col_name = "time",
            permute = FALSE
        ),
        "Run GCalignR"
    )

    quiet_result <- NULL
    options_before_quiet <- pbapply::pboptions()
    expect_output(
        quiet_result <- align_chromatograms(
            peak_subset,
            rt_col_name = "time",
            permute = FALSE,
            verbose = FALSE
        ),
        NA
    )

    expect_identical(pbapply::pboptions(), options_before_quiet)
    expect_equal(quiet_result$aligned_list, verbose_result$aligned_list)
})

test_that("quiet reference selection restores pbapply options", {
    data("peak_data", package = "GCalignR")
    peak_subset <- lapply(peak_data[1:3], function(x) x[1:50, ])
    old_options <- pbapply::pboptions(type = "txt")
    on.exit(pbapply::pboptions(old_options))
    options_before <- pbapply::pboptions()

    expect_output(
        choose_optimal_reference(peak_subset, rt_col_name = "time", verbose = FALSE),
        NA
    )
    expect_identical(pbapply::pboptions(), options_before)
})
