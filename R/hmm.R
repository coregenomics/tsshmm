#' @import GenomicRanges
#' @importFrom BiocGenerics relist start score
#' @importFrom S4Vectors endoapply mendoapply mcols `mcols<-` runValue
#'     queryHits subjectHits to
#' @importFrom methods as is
#' @importFrom stats aggregate
NULL

encode <- function(signal, bg, windows) {
    stopifnot(is(signal, "GRanges"))
    stopifnot(is(bg, "GRanges"))
    stopifnot(is(windows, "CompressedGRangesList"))

    win_signal <- scoreOverlapWindows(windows, signal)
    win_bg <- scoreOverlapWindows(windows, bg)

    ## Our observation states are:
    ## 0) No signal state (signal == 0)
    ## 1) Enriched        (signal > bg)
    ## 2) Depleted  (bg >= signal > 0)
    ##
    no_sig <- win_signal == 0
    enriched <- win_signal > win_bg
    depleted <- win_bg >= win_signal

    obs <- relist(
        rep(NA, sum(lengths(win_signal))),
        win_signal)
    obs <- as(obs, "IntegerList")
    obs[depleted] <- 2L
    obs[enriched] <- 1L
    obs[no_sig] <- 0L
    stopifnot(sum(unlist(is.na(obs))) == 0)

    obs
}

check_valid_hmm_reads <- function(gr) {
    stopifnot(is(gr, "GRanges"))
    if (length(gr) == 0L)
        return()
    stopifnot(any(strand(gr) != "*"))
    stopifnot(! is.null(score(gr)))
    stopifnot(all(score(gr) >= 0L))
}

stranded <- function(x) split(x, strand(x))[-3] # -3 removes "*" strand.

hmm <- function(signal, bg, ranges) {
    check_valid_hmm_reads(signal)
    check_valid_hmm_reads(bg)
    stopifnot(is(ranges, "GRanges"))
    mendoapply(hmm_by_strand, stranded(signal), stranded(bg),
               MoreArgs = list(ranges = ranges))
}

## Rely on annotations for now.  Based on the data below, more careful work is
## required to adjust the methodology to be annotation free.  This is because
## choosing thresholds for how far apart GRO-cap counts can be to consider them
## as part of a promoter is tricky.  Most likely, one would need to use the HMM
## calls to inform those thresholds.
##
## ENCODE inter-enhancer distances on chr22:
##
## > mcols(distanceToNearest(ranges))$distance %>%
##   quantile(seq(0.05, 1, 0.05)) %>% as("integer")
##
##  5%   10%   15%   20%   25%   30%   35%   40%   45%   50%
##   2     5    10    21    35    54    74    97   126   160
## 55%   60%   65%   70%   75%   80%   85%   90%   95%  100%
## 195   238   291   351   427   536   714  1036  1765 27816
##
## Looking at + strand counts subset to those enhancers:
##
## > dist <- distanceToNearest(signal)
## > mcols(dist)$distance %>% quantile(seq(0.05, 1, 0.05)) %>% as("integer")
##
##  5%    10%    15%    20%    25%    30%    35%    40%    45%    50%
##   0      0      0      0      0      0      1      1      1      2
## 55%    60%    65%    70%    75%    80%    85%    90%    95%   100%
##   3      4      6      9     15     25     47    122    641 191170
hmm_by_strand <- function(signal, bg, ranges) {
    strand <- runValue(strand(signal))
    if (length(strand) == 0L)
        return(GRanges())
    stopifnot(length(strand) == 1L)
    stopifnot(all(strand == runValue(strand(bg))))
    ## No signal at all in 80% of the annotated chr22+ enhancers regions!
    ## Therefore save calculation time by removing ranges without any counts.
    ##
    ## > (scoreOverlaps(ranges, signal) > 0) %>% table() %>%
    ##   `/`(length(ranges) / 100) %>% print(digits = 3)
    ##
    ## FALSE  TRUE
    ##  80.3  19.7
    ##
    has_signal <- scoreOverlaps(ranges, signal) > 0
    if (any(has_signal == 0))
        message(sprintf("Dropping %.1f%% of regions with no signal.",
                        100 * (1 - sum(has_signal) / length(ranges))))
    windows <- tile(ranges[has_signal], width = 10)
    if (strand == "-")
        windows <- endoapply(windows, rev)
    observations <- encode(signal, bg, windows)
    is_promoter <- endoapply(observations, viterbi) > 0

    gr <- reduce(unlist(windows)[unlist(is_promoter)])
    strand(gr) <- strand
    gr
}

replace_unstranded <- function (gr) {
    idx <- strand(gr) == "*"
    if (length(idx) == 0L)
        return(gr)
    sort(c(
        gr[! idx],
        `strand<-`(gr[idx], value = "+"),
        `strand<-`(gr[idx], value = "-")))
}

tss <- function(signal, ranges) {
    ## Find the highest peak within each region using 2 bp windows.  Tiling 2
    ## bp windows is inefficient, so instead use a rolling max in the C layer.
    stopifnot(is(ranges, "GRanges"))
    stopifnot(is(signal, "GRanges"))
    if (any(strand(signal) == "*"))
        stop("signal must be stranded")
    ranges <- replace_unstranded(ranges)
    ol <- findOverlaps(ranges, signal)
    groups <- queryHits(ol)
    indices_signal <- subjectHits(ol)
    starts_signal <-
        start(signal[subjectHits(ol)]) -
        start(ranges[queryHits(ol)])
    scores_signal <- score(signal[subjectHits(ol)])
    ## Down convert numeric to integer.
    if (! is(scores_signal, "integer"))
        scores_signal <- as.integer(scores_signal)
    n_groups <- length(unique(queryHits(ol)))
    indices_peak <- vector("integer", n_groups)
    scores_peak <- vector("integer", n_groups)
    .Call(
        C_tss, PACKAGE = "tsshmm", indices_peak, groups, indices_signal,
        starts_signal, scores_signal)
    signal[indices_peak]
}

scoreOverlaps <- function(gr, reads) {
    ol <- findOverlaps(gr, reads)
    if (length(ol) == 0L)
        return(vector("integer", length(gr)))
    mcols(ol)$score <- score(reads[to(ol)])
    df <- as(ol, "DataFrame")
    df <- aggregate(score ~ queryHits, df, max)
    score <- vector("integer", length(gr))
    score[df$queryHits] <- df$score
    score
}

scoreOverlapWindows <- function(grl, reads)
    relist(scoreOverlaps(unlist(grl), reads), grl)
