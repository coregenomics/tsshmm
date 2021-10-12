tol <- 10 ^ (-1 * 3:8)
## These expected values are calculated in the vignette appendix.
x <- c(37, 48, 59, 69, 80, 91)
stopifnot(length(tol) == length(x))
n_low <- floor(x / 10) * 10
n_low <- ifelse(n_low == x, floor((x - 1) / 10) * 10, n_low)
stopifnot(all(n_low < x))
n_high <- n_low + 10
n_high <- ifelse(n_high == x, n_high + 10, n_high)
stopifnot(all(n_high > x))

context("viterbi helpers")

test_that("prom_dist is valid for several tolerances", {
    model <- TSSHMM()
    for (i in seq_along(tol)) {
        expect_equal(prom_dist(model, tol = tol[i]), x[i])
    }
})

test_that("prom_dist detects tolerance too stringent", {
    model <- TSSHMM()
    for (i in seq_along(tol)) {
        expect_error(prom_dist(model, tol = tol[i], n = n_low[i]))
        expect_silent(prom_dist(model, tol = tol[i], n = n_high[i]))
    }
})

test_that("stranded always splits into + and - levels", {
    expect_equal(lengths(stranded(GRanges("chr1:100:+"))), c("+" = 1, "-" = 0))
    expect_equal(lengths(stranded(GRanges("chr1:100:-"))), c("+" = 0, "-" = 1))
    expect_equal(lengths(stranded(
        GRanges(c("chr1:100:+", "chr1:200:-")))),          c("+" = 1, "-" = 1))
    expect_equal(lengths(stranded(GRanges("chr1:100:*"))), c("+" = 0, "-" = 0))
    expect_equal(lengths(stranded(GRanges())),             c("+" = 0, "-" = 0))
})

context("viterbi")

test_that("viterbi_by_strand joins signal regions", {
    model <- TSSHMM()
    promoters <- viterbi_by_strand(model, signal, bg, ranges)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:101-300:+",
                           "chr2:201-260:+",
                           "chr2:281-450:+")))
    gr_mask <- unlist(bases)[as.logical(mask)]
    n_ol <- sum(countOverlaps(gr_mask, promoters))
    expect_gt(n_ol / sum(mask), 0.8)
})

test_that("viterbi_by_strand notifies about missing regions", {
    model <- TSSHMM()
    expect_silent(viterbi_by_strand(model, signal, bg, ranges))
    expect_message(viterbi_by_strand(model, signal, bg, c(ranges, empty)),
                   "Dropping 33[.]3%")
})

test_that("viterbi runs on both strands", {
    model <- TSSHMM()
    promoters <- viterbi(model, signal, bg)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:102-300:+",
                           "chr2:202-260:+",
                           "chr2:281-449:+")))
    strand(signal) <- "-"
    promoters <- viterbi(model, signal, bg)
    promoters_no_metadata <- promoters
    mcols(promoters_no_metadata) <- NULL
    expect_equal(promoters_no_metadata,
                 GRanges(c("chr1:96-283:-",
                           "chr1:285-293:-",
                           "chr2:219-437:-",
                           "chr2:439-447:-")))
})
