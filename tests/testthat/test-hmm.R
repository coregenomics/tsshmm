ranges <- GRanges(c("chr1:101-300",
                    "chr2:201-450"))
n <- sum(width(ranges))
bases <- tile(ranges, width = 1)
strand(bases) <- "+"
set.seed(123)
mask <- rbinom(n, 1, 0.3)
score_sig <- rnbinom(n, 1, 0.5) * mask
setScore <- function(gr, score) {
    scored <- `score<-`(unlist(gr), value = score)
    subset(scored, score(scored) > 0)
}
signal <- setScore(bases, score_sig)
score_bg <- rnbinom(n, 0.7, 0.5)
bg <- setScore(bases, score_bg)
empty <- GRanges("chr1:600-650")

context("hmm")

test_that("hmm joins signal regions", {
    promoters <- hmm(signal, bg, ranges)
    expect_equal(promoters, GRanges(c("chr1:111-300",
                                      "chr2:211-260",
                                      "chr2:281-450")))
    gr_mask <- unlist(bases)[as.logical(mask)]
    n_ol <- sum(countOverlaps(gr_mask, promoters))
    expect_gt(n_ol / sum(mask), 0.8)
})

test_that("hmm notifies about missing regions", {
    expect_silent(hmm(signal, bg, ranges))
    expect_message(hmm(signal, bg, c(ranges, empty)), "Dropping 33[.]3%")
})

context("encode")

test_that("encode checks inputs", {
    windows <- tile(ranges, width = 10)
    expect_silent(encode(signal, bg, windows))
    expect_error(encode(NULL, bg, windows), "signal")
    expect_error(encode(signal, NULL, windows), "bg")
    expect_error(encode(signal, bg, NULL), "windows")
})

context("scoreOverlaps")

test_that("scoreOverlaps min-maxes positive and negative counts", {
    gr <- c(ranges, empty)
    ol <- findOverlaps(gr, signal)
    df <- data.frame(i = queryHits(ol),
                     score = score(signal[subjectHits(ol)]))
    df <- aggregate(score ~ ., data = df, max)
    expect_equal(scoreOverlaps(gr, signal), c(df$score, 0))
    score(signal) <- -score(signal)
    expect_equal(scoreOverlaps(gr, signal), -c(df$score, 0))
})

context("tss")

test_that("tss returns peak values and positions", {
    gr <- c(ranges, empty)
    ol <- findOverlaps(gr, signal)
    df <- data.frame(i = queryHits(ol),
                     score = score(signal[subjectHits(ol)]))
    df <- aggregate(score ~ ., data = df, max)
    tss <- tss(signal, gr)
    expect_equal(score(tss), df$score)
})

test_that("tss handles non-integer scores", {
    signal <- GRanges(c("chr1:100",
                        "chr1:200"),
                      score = c(9L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
    score(signal) <- c(9.0, 9.0)
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
})

test_that("tss breaks ties using look head and look behind", {
    ## Should choose first peak encountered when there is no adjacency.
    signal <- GRanges(c("chr1:100",
                        "chr1:200"),
                      score = c(9L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
    ## Look head.
    signal <- GRanges(c("chr1:100",
                        "chr1:200",
                        "chr1:201"),
                      score = c(9L, 9L, 1L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[2])
    ## Look behind.
    signal <- GRanges(c("chr1:100",
                        "chr1:200",
                        "chr1:201"),
                      score = c(9L, 1L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[3])
})
