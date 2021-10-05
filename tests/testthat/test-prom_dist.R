context("viterbi helpers")

tol <- 10 ^ (-1 * 3:8)
## These expected values are calculated in the vignette appendix.
x <- c(28, 39, 50, 60, 71, 82)
stopifnot(length(tol) == length(x))
n_low <- floor(x / 10) * 10
n_low <- ifelse(n_low == x, floor((x - 1) / 10) * 10, n_low)
stopifnot(all(n_low < x))
n_high <- n_low + 10
n_high <- ifelse(n_high == x, n_high + 10, n_high)
stopifnot(all(n_high > x))

test_that("prom_dist is valid for several tolerances", {
    model <- new("TSSHMM")
    for (i in seq_along(tol)) {
        expect_equal(prom_dist(model, tol = tol[i]), x[i])
    }
})

test_that("prom_dist detects tolerance too stringent", {
    model <- new("TSSHMM")
    for (i in seq_along(tol)) {
        expect_error(prom_dist(model, tol = tol[i], n = n_low[i]))
        expect_silent(prom_dist(model, tol = tol[i], n = n_high[i]))
    }
})
