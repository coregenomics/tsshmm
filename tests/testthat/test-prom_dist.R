context("viterbi helpers")

test_that("prom_dist is valid for several tolerances", {
    model <- new("TSSHMM")
    ## These expected values are calculated in the vignette appendix.
    expect_equal(prom_dist(model, tol = 10^-3), 28)
    expect_equal(prom_dist(model, tol = 10^-4), 39)
    expect_equal(prom_dist(model, tol = 10^-5), 50)
    expect_equal(prom_dist(model, tol = 10^-6), 60)
    expect_equal(prom_dist(model, tol = 10^-7), 71)
    expect_equal(prom_dist(model, tol = 10^-8), 82)    
})
