context("parameters")

test_that("TSSHMM states and emissions can be read with dim", {
    model <- TSSHMM(bg_genebody = FALSE)
    expect_equal(dim(model), c(states = 7, emissions = 3))
    model <- TSSHMM(bg_genebody = TRUE)
    expect_equal(dim(model), c(states = 8, emissions = 3))
})

test_that("TSSHMM parameters can be saved and loaded", {
    model <- TSSHMM()
    params_orig <- parameters(model)
    params_new <- params_orig
    params_new$trans["B", ] <- c(0.9, 0.05, NA, NA, 0.05, NA, NA)
    stopifnot(all(rowSums(params_new$trans, na.rm = TRUE) == 1))
    params_new$emis[paste0("P", 1:3), ] <- rep(c(0.2, 0.4, 0.4), each = 3)
    stopifnot(all(rowSums(params_new$emis) == 1))
    parameters(model) <- params_new
    expect_equal(parameters(model), params_new)
})

test_that("params_names generates column-major concatenated dimnames", {
    mat <- matrix(nrow = 10, ncol = 5,
                  dimnames = list(LETTERS[1:10], letters[11:15]))
    expect_equal(params_names(mat, 1:5), paste0(LETTERS[1:5], ".k"))
    expect_equal(params_names(mat, 11:15), paste0(LETTERS[1:5], ".l"))
    expect_equal(params_names(mat, 6:15), c(paste0(LETTERS[6:10], ".k"),
                                            paste0(LETTERS[1:5], ".l")))
    expect_equal(params_names(mat, 26:30), paste0(LETTERS[6:10], ".m"))
    expect_equal(params_names(mat, (length(mat)-4):length(mat)),
                 paste0(LETTERS[6:10], ".o"))
    expect_equal(params_names(mat, (length(mat)+1):(length(mat)+5)),
                 rep("NA.NA", 5))
})

test_that("params_idx", {
    # Regression test of externalptr of the second model bizarrely
    # being shared with the first model.  This test was fixed by no
    # longer using externalptr anywhere in the code and moving all the
    # low-level C mode logic into R.
    model <- TSSHMM(bg_genebody = FALSE)
    n_params_procap <- length(unlist(params_idx(model)))
    model <- TSSHMM(bg_genebody = TRUE)
    n_params_proseq <- length(unlist(params_idx(model)))
    expect_gt(n_params_proseq,
              n_params_procap)
})
