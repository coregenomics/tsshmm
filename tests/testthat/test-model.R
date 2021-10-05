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

test_that("TSSHMM parameters can be saved and loaded", {
    model <- new("TSSHMM")
    params_orig <- parameters(model)
    params_new <- params_orig
    params_new$trans["B", ] <- c(0.9, 0.05, NA, NA, 0.05, NA, NA)
    stopifnot(all(rowSums(params_new$trans, na.rm = TRUE) == 1))
    params_new$emis[paste0("P", 1:3), ] <- rep(c(0.2, 0.4, 0.4), each = 3)
    stopifnot(all(rowSums(params_new$emis) == 1))
    parameters(model) <- params_new
    expect_equal(parameters(model), params_new)
})
