context("scoreOverlaps")

test_that("scoreOverlaps finds max score per range", {
    gr <- c(ranges, empty)
    ol <- findOverlaps(gr, signal)
    df <- data.frame(i = queryHits(ol),
                     score = score(signal[subjectHits(ol)]))
    df <- aggregate(score ~ ., data = df, max)
    expect_equal(scoreOverlaps(gr, signal), c(df$score, 0))
})

test_that("scoreOverlaps handles empty input", {
    expect_equal(scoreOverlaps(c(empty, empty), signal), c(0L, 0L))
})
