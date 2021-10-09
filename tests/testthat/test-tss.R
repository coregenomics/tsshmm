context("tss helpers")

test_that("replace_unstranded ignores stranded regions", {
    gr <- GRanges(c("chr:100:+", "chr:200:-"))
    expect_equal(replace_unstranded(gr), gr)
    gr <- GRanges(c("chr:100:-", "chr:200:+"))
    expect_equal(replace_unstranded(gr), sort(gr))
})

test_that("replace_unstranded splits unstranded into + and -", {
    expect_equal(replace_unstranded(GRanges("chr:100")),
                 GRanges(c("chr:100:+", "chr:100:-")))
    expect_equal(replace_unstranded(GRanges(c("chr:100", "chr:200:+"))),
                 sort(GRanges(c("chr:100:+", "chr:100:-", "chr:200:+"))))
})

test_that("replace_unstranded ignores empty GRanges", {
    gr <- GRanges()
    expect_equal(replace_unstranded(gr), gr)
})

context("tss")

test_that("tss checks inputs", {
    gr <- GRanges("chr:51-150")
    expect_error(tss(GRanges("chr:100"), gr), "strand")
    expect_silent(tss(GRanges("chr:100:+"), gr))
    expect_silent(tss(GRanges("chr:100:-"), gr))
    expect_error(tss(GRanges("chr:100:-"), NULL), "GRanges")
    expect_error(tss(NULL, gr), "GRanges")
})

test_that("tss requires sorted signal", {
    gr <- GRanges("chr:51-150")
    signal <- GRanges(c("chr:150:+", "chr:100:+"))
    expect_error(tss(signal, gr), "isSorted")
    expect_silent(tss(sort(signal), gr))
})

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
    signal <- GRanges(c("chr1:100:+",
                        "chr1:200:+"),
                      score = c(9L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
    score(signal) <- c(9.0, 9.0)
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
})

test_that("tss breaks ties using look head and look behind", {
    ## Should choose first peak encountered when there is no adjacency.
    signal <- GRanges(c("chr1:100:+",
                        "chr1:200:+"),
                      score = c(9L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[1])
    ## Look ahead.
    signal <- GRanges(c("chr1:100:+",
                        "chr1:200:+",
                        "chr1:201:+"),
                      score = c(9L, 9L, 1L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[2])
    ## Look behind.
    signal <- GRanges(c("chr1:100:+",
                        "chr1:200:+",
                        "chr1:201:+"),
                      score = c(9L, 1L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[3])
    ## Look ahead and look behind.
    signal <- GRanges(c("chr1:100:+",
                        "chr1:101:+",
                        "chr1:199:+",
                        "chr1:200:+",
                        "chr1:201:+"),
                      score = c(2L, 9L, 2L, 9L, 1L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[4])
})

test_that("tss breaks ties strand-specifically", {
    ## Choose the first peak on the positive strand, and last peak on
    ## the negative strand.
    signal <- GRanges(c(
            "chr1:100:+", ## first peak on the positive strand that we want.
            "chr1:199:+", ## last peak on the positive strand.
            "chr1:101:-", ## first peak on the negative strand.
            "chr1:200:-"  ## last peak on the negative strand that we want.
    ), score = c(9L, 9L, 9L, 9L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, signal[c(1, 4)])
})

test_that("tss returns strand-specific peaks", {
    signal <- GRanges(c("chr1:100:+",
                        "chr1:200:+",
                        "chr1:150:-",
                        "chr1:250:-"),
                      score = c(1L, 9L, 9L, 1L))
    tss <- tss(signal, range(signal))
    expect_equal(tss, sort(signal[2:3]))
    tss <- tss(signal, range(signal, ignore.strand = TRUE))
    expect_equal(tss, sort(signal[2:3]))
})

test_that("tss optionally returns Pairs with query range", {
    gr <- c(ranges, empty)
    ol <- findOverlaps(gr, signal)
    df <- data.frame(i = queryHits(ol),
                     score = score(signal[subjectHits(ol)]))
    df <- aggregate(score ~ ., data = df, max)
    tss <- tss(signal, gr)
    pairs <- tss(signal, gr, pairs = TRUE)
    expect_type(pairs, "S4")
    expect_s4_class(pairs, "Pairs")
    expect_equal(S4Vectors::first(pairs), tss)
    expect_equal(second(pairs), `strand<-`(ranges, value = "+"))
})
