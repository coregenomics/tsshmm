futile.logger::flog.threshold(futile.logger::ERROR)
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

context("hmm helpers")

test_that("check_valid_hmm_reads flags invalid input", {
    expect_error(check_valid_hmm_reads(NULL), "GRanges")
    expect_error(check_valid_hmm_reads(GRanges("chr1:100:*")), "strand.+[*]")
    expect_error(check_valid_hmm_reads(GRanges("chr1:100:+")), "score")
    expect_error(check_valid_hmm_reads(GRanges("chr1:100:+", score = -10)),
                 "score.*>=.*0")
    expect_silent(check_valid_hmm_reads(GRanges("chr1:100:+", score = 0)))
    expect_silent(check_valid_hmm_reads(GRanges("chr1:100:-", score = 1.5)))
    expect_silent(check_valid_hmm_reads(GRanges()))
})

test_that("stranded always splits into + and - levels", {
    expect_equal(lengths(stranded(GRanges("chr1:100:+"))), c("+" = 1, "-" = 0))
    expect_equal(lengths(stranded(GRanges("chr1:100:-"))), c("+" = 0, "-" = 1))
    expect_equal(lengths(stranded(
        GRanges(c("chr1:100:+", "chr1:200:-")))),          c("+" = 1, "-" = 1))
    expect_equal(lengths(stranded(GRanges("chr1:100:*"))), c("+" = 0, "-" = 0))
    expect_equal(lengths(stranded(GRanges())),             c("+" = 0, "-" = 0))
})

test_that("tile_with_rev is equivalent to endoapply(tile(...), rev)", {
    tiled <- tile(ranges, width = 10)
    expect_equal(tile_with_rev(ranges, 10, rev = FALSE), tiled)
    expect_equal(tile_with_rev(ranges, 10, rev = TRUE), endoapply(tiled, rev))
})

context("hmm")

test_that("viterbi_by_strand joins signal regions", {
    model <- new("TSSHMM")
    promoters <- viterbi_by_strand(model, signal, bg, ranges)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:111-300:+",
                           "chr2:231-260:+",
                           "chr2:281-410:+",
                           "chr2:431-450:+")))
    gr_mask <- unlist(bases)[as.logical(mask)]
    n_ol <- sum(countOverlaps(gr_mask, promoters))
    expect_gt(n_ol / sum(mask), 0.8)
})

test_that("viterbi_by_strand notifies about missing regions", {
    model <- new("TSSHMM")
    expect_silent(viterbi_by_strand(model, signal, bg, ranges))
    expect_message(viterbi_by_strand(model, signal, bg, c(ranges, empty)),
                   "Dropping 33[.]3%")
})

test_that("viterbi runs on both strands", {
    model <- new("TSSHMM")
    promoters <- viterbi(model, signal, bg)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:111-250:+",
                           "chr1:261-290:+",
                           "chr2:231-260:+",
                           "chr2:281-409:+",
                           "chr2:420-449:+")))
    strand(signal) <- "-"
    promoters <- viterbi(model, signal, bg)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:97-106:-",
                           "chr1:108-135:-",
                           "chr1:146-165:-",
                           "chr1:167-214:-",
                           "chr1:226-283:-",
                           "chr2:219-268:-",
                           "chr2:279-337:-",
                           "chr2:348-437:-")))
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
    expect_equal(first(pairs), tss)
    expect_equal(second(pairs), `strand<-`(ranges, value = "+"))
})
