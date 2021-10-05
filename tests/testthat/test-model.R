context("model")

test_that("TSSHMM model initializes with bg_proseq toggle", {
    model_procap <- new("TSSHMM", bg_proseq = FALSE)
    expect_equal(dim(parameters(model_procap)$emis), c(7, 3))
    model_proseq <- new("TSSHMM", bg_proseq = TRUE)
    expect_equal(dim(parameters(model_proseq)$emis), c(8, 3))
})

test_that("TSSHMM model flags invalid bg_proseq input", {
    expect_error(new("TSSHMM", bg_proseq = 1.0))
    expect_error(new("TSSHMM", bg_proseq = -1))
    expect_error(new("TSSHMM", bg_proseq = 0))
    expect_error(new("TSSHMM", bg_proseq = "yes"))
    expect_error(new("TSSHMM", bg_proseq = "no"))
    expect_silent(new("TSSHMM", bg_proseq = TRUE))
    expect_silent(new("TSSHMM", bg_proseq = FALSE))
})
