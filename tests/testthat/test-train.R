context("train")

test_that("train does not change model for empty data", {
    model <- new("TSSHMM")
    params <- parameters(model)
    batches <- list()
    updates <- train(model, batches)
    expect_equal(parameters(model), params)
    expect_equal(nrow(updates), 1)
})

test_that("train converges to true parameters", {
    model <- new("TSSHMM")
    params_actual <- parameters(model)

    params_start <- params_actual
    params_start$trans["B", c("B", "N1", "P1")] = c(0.9, 0.05, 0.05)
    stopifnot(all(rowSums(params_start$trans, na.rm = TRUE) == 1))

    obs <- matrix(-1L, nrow = 100, ncol = 1e4)
    set.seed(123)
    .Call(C_simulate, PACKAGE = "tsshmm", obs, model@external_pointer)
    ## https://stackoverflow.com/a/6821395
    list_from_matrix <- function(x) lapply(seq_len(ncol(x)), function(i) x[,i])
    batches <- list(IntegerList(list_from_matrix(t(obs))))
    updates <- train(model, batches)
    diffs <-
        function(params) abs(params$trans["B", c("B", "N1", "P1")] -
                             parameters(model)$trans["B", c("B", "N1", "P1")])
    expect_true(all(diffs(params_actual) * 10 < diffs(params_start)))
})
