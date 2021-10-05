context("validate reads")

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
