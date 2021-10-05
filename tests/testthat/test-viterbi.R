context("viterbi helpers")

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
