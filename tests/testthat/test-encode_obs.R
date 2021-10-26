context("encode_obs helpers")

test_that("tile_with_rev is equivalent to endoapply(tile(...), rev)", {
    tiled <- tile(ranges, width = 10)
    expect_equal(tile_with_rev(ranges, 10, rev = FALSE), tiled)
    expect_equal(tile_with_rev(ranges, 10, rev = TRUE), endoapply(tiled, rev))
})

context("encode_obs")

test_that("encode_obs allows empty input", {
    obs <- encode_obs(GRanges(), GRanges())
    expect_equal(obs, RleList())
})

test_that("encode_obs ignores rows with only zeros", {
    obs <- IntegerList()
    n_cols <- 6
    n_rows <- 10
    signal <- unlist(tile(
        GRanges(paste0("chr", 1:n_rows, ":100-1000:+")), n_cols))
    score(signal) <- seq_along(signal)
    obs <- encode_obs(signal, GRanges(), nrow = n_cols)
    expect_equal(length(obs), n_rows)
    n_rows_zeros <- 3
    only_zeros <- sample(n_rows, n_rows_zeros)
    mask <- rep(1:n_rows, each = n_cols) %in% only_zeros
    score(signal[mask]) <- 0
    obs <- encode_obs(signal, GRanges(), nrow = n_cols)
    expect_equal(length(obs), n_rows - n_rows_zeros)
})

test_that("encode_obs segments up to 1e4 width 10 windows", {
    signal <- GRanges("chr:1-1e6:+", score = 2)
    obs <- encode_obs(signal, GRanges())
    expect_equal(length(obs), 10)
    expect_equal(elementNROWS(obs), rep(1e4, 10))
})

test_that("encode_obs with sample shuffles by rows", {
    gr <- tile(GRanges("chr:1-1e6:+"), 1e4)[[1]]
    set.seed(123)
    scores <- function(theta)
        `score<-`(gr, value = MASS::rnegbin(length(gr), theta = theta))
    signal <- scores(0.1)
    bg <- scores(1)
    obs <- list()
    seeds <- c(123, 456, 789)
    for (i in seq_along(seeds)) {
        obs[[i]] <- encode_obs(signal, bg)
        set.seed(seeds[i])
        obs[[i]] <- sample(obs[[i]])
    }
    digests <- function(x) sapply(x, digest::digest, USE.NAMES = FALSE)
    for (i in 1:(length(seeds) - 1)) {
        sort_batch <- function(x) apply(t(as.matrix(x)), 1, sort)
        expect_equal(sort_batch(obs[[i]]),
                     sort_batch(obs[[i+1]]))
        expect_false(digest::digest(obs[[i]]) ==
                     digest::digest(obs[[i+1]]))
    }
})
