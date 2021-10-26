context("train")

test_that("train does not change model for empty data", {
    model_orig <- TSSHMM()
    params <- parameters(model_orig)
    obs <- RleList()
    model <- train(model_orig, obs)
    expect_equal(parameters(model), params)
    expect_identical(model, model_orig)
})

test_that("train converges to true parameters", {
    model <- TSSHMM()
    params_actual <- parameters(model)

    params_start <- params_actual
    params_start$trans["B", c("B", "N1", "P1")] = c(0.9, 0.05, 0.05)
    stopifnot(all(rowSums(params_start$trans, na.rm = TRUE) == 1))

    set.seed(123)
    obs <-
        .Call(C_simulate, PACKAGE = "tsshmm", c(100L, 1e4L),
              transitions(model), emissions(model),
              emissions_tied(model), start(model))
    expect_equal(digest::digest(obs), "f1681b80ab2b8ad0512a544cebfa4db0")
    ## https://stackoverflow.com/a/6821395
    list_from_matrix <- function(x) lapply(seq_len(ncol(x)), function(i) x[,i])
    obs <- as(list_from_matrix(t(obs)), "RleList")
    parameters(model) <- params_start
    set.seed(456)
    model <- train(model, obs)
    diffs <-
        function(params) abs(params$trans["B", c("B", "N1", "P1")] -
                             parameters(model)$trans["B", c("B", "N1", "P1")])
    expect_true(all(diffs(params_actual) * 10 < diffs(params_start)))
})
