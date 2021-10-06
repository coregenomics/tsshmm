context("create_batches helpers")

test_that("tile_with_rev is equivalent to endoapply(tile(...), rev)", {
    tiled <- tile(ranges, width = 10)
    expect_equal(tile_with_rev(ranges, 10, rev = FALSE), tiled)
    expect_equal(tile_with_rev(ranges, 10, rev = TRUE), endoapply(tiled, rev))
})

context("create_batches")

test_that("create_batches allows empty input", {
    batches <- create_batches(GRanges(), GRanges())
    expect_equal(batches, list())
})

test_that("create_batches segments up to 1e4 width 10 windows", {
    signal <- GRanges("chr:1-1e6:+", score = 2)
    batches <- create_batches(signal, GRanges())
    expect_equal(elementNROWS(batches), 10)
    expect_equal(elementNROWS(batches[[1]]),
                 rep(1e4, 10))
    batches <- create_batches(signal, GRanges(), nrow = 3)
    elementNROWS <- elementNROWS(batches)
    expect_equal(sum(elementNROWS), 10)
    expect_equal(length(elementNROWS), 4)
})

test_that("create_batches shuffles by rows and not by batch", {
    gr <- tile(GRanges("chr:1-1e6:+"), 1e4)[[1]]
    set.seed(123)
    scores <- function(theta)
        `score<-`(gr, value = MASS::rnegbin(length(gr), theta = theta))
    signal <- scores(0.1)
    bg <- scores(1)
    batches <- list()
    seeds <- c(123, 456, 789)
    for (i in seq_along(seeds)) {
        batches[[i]] <- create_batches(signal, bg, seed = seeds[i], nrow = 3)
    }
    digests <- function(x) sapply(x, digest::digest, USE.NAMES = FALSE)
    for (i in 1:(length(seeds) - 1)) {
        sort_batch <- function(x) apply(t(as.matrix(unlist(List(x)))), 1, sort)
        expect_equal(sort_batch(batches[[i]]),
                     sort_batch(batches[[i+1]]))
        n_diffs <- length(setdiff(digests(batches[[i]]),
                                  digests(batches[[i+1]])))
        expect_equal(n_diffs, length(batches[[i]]))
    }
})
