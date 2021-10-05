#' @import GenomicRanges
#' @importFrom BiocGenerics relist start score
#' @importFrom IRanges IRanges PartitioningByWidth width
#' @importFrom S4Vectors elementNROWS endoapply isSorted mendoapply mcols
#'     "mcols<-" Pairs runValue queryHits subjectHits to
#' @importFrom futile.logger flog.debug
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
viterbi_by_strand <- function(model, signal, bg, ranges) {
    flog.debug("Entered viterbi_by_strand()")
    flog.debug("Checking inputs")
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
    flog.debug("Trimming regions with no signal")
    has_signal <- scoreOverlaps(ranges, signal) > 0
    if (any(has_signal == 0))
        message(sprintf("Dropping %.1f%% of regions with no signal.",
                        100 * (1 - sum(has_signal) / length(ranges))))
    flog.debug("Generating windows")
    windows <- tile_with_rev(ranges[has_signal], 10, rev = (strand == "-"))
    flog.debug("Encoding observations as enriched, depleted, or background")
    observations <- encode(signal, bg, windows)
    flog.debug("Running Viterbi")
    states <- viterbi_low_level(model, observations)
    flog.debug("Reducing discovered regions")
    is_promoter <- states > 0

    gr <- unlist(windows)
    strand(gr) <- strand
    gr$states <- unlist(states)
    ## Several hidden states may compose a promoter region, therefore preserve
    ## all states using the with.revmap option.
    gr_reduced <- reduce(gr[unlist(is_promoter)], with.revmap = TRUE)
    flog.debug("Preserving hidden states")
    gr_reduced$states <- relist(
        gr$states[unlist(is_promoter)][unlist(gr_reduced$revmap)],
        gr_reduced$revmap)
    gr_reduced$revmap <- NULL
    flog.debug("Exiting viterbi_by_strand()")
    gr_reduced
}

## Running `grl <- tile(gr)` followed by `endoapply(grl, rev)` is extremely
## slow (tile() is a few seconds, but rev() is about an hour!), therefore we
## need to support this narrow case for efficient reversed tiles by making
## slight modifications to the original GenomicRanges::tile() method and it's
## inner IRanges::tile() call.
##
## The few lines changed from this methods have been noted below.  Extraneous
## lines from the original functions have been removed because don't apply to
## our narrow use case: in this package, we only care about the `width`
## argument, and not the `n` argument.
##
## Perhaps a patch to GenomicRanges::tile(..., reverse = FALSE) with this
## efficient feature can be contributed upstream because it would be generally
## useful for any nascent RNA or DNA sequencing analysis.  If and when such a
## patch is accepted, then this we can go back to using GenomicRanges::tile()
## and drop this workaround.
tile_with_rev <- function(x, width, rev) {
    ## Begin GenomicRanges::tile() [GenomicRanges version 1.45.0]
    seqnames <- seqnames(x)
    strand <- strand(x)
    x <- ranges(x)
    ## Begin IRanges::tile() [IRanges version 2.27.0]
    n <- ceiling(width(x) / width)
    width <- IRanges::width(x) / n
    ## Begin refactored lines to simplify our change.
    ir <- IRanges(rep(1L, length(n)), width = n)
    vec.ref <- sequence(width(ir), from = start(ir))
    if (rev) {
        ## Begin lines we're functionally changing.
        vec.end <- sequence(width(ir), from = end(ir), by = -1L)
        ## End lines we're functionally changing.
    } else {
        vec.end <- vec.ref
    }
    tile.end <- floor(vec.end * rep(width, n))
    tile.ref <- floor(vec.ref * rep(width, n))
    tile.end.abs <- tile.end + rep(start(x), n) - 1L
    tile.width <- S4Vectors:::diffWithInitialZero(as.integer(tile.ref))
    ## End refactored lines to simplify our change.
    p <- IRanges:::PartitioningByWidth(n, names = names(x))
    tile.width[start(p)] <- tile.ref[start(p)]
    tiles <- relist(IRanges(width = tile.width,
                            end = tile.end.abs), p)
    ## End: IRanges::tile()
    gr <- GRanges(rep(seqnames, elementNROWS(tiles)), unlist(tiles),
                  rep(strand, elementNROWS(tiles)))
    relist(gr, tiles)
    ## End GenomicRanges::tile()
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

#' Find promoter peaks using 3 basepair tie-breaking.
#'
#' The original method for extracting the TSS was using 2 basepair windows, but
#' note that this detail was not documented in
#' \insertCite{core_analysis_2014}{tsshmm}.  Preserving the intention of that
#' method to use neighboring counts to break ties, this method uses 3 basepair
#' overlapping windows, where the preceding and succeeding basepair counts act
#' as a "bonus" value to help break ties.
#'
#' @param signal Stranded, single base \code{GRanges} with integer score.
#' @param ranges \code{GRanges} to limit signal search.
#' @param pairs Whether to embed the \code{GRanges} into a \code{Pairs} object
#'     with the stranded range in which the TSS was found.
#' @return Stranded, single base \code{GRanges} with integer score.  If
#' @examples
#' # When counts are equal, TSS returns first maximum.
#' signal <- GRanges(c("chr1:100:+", "chr1:200:+"), score = c(9L, 9L))
#' signal
#' tss <- tss(signal, range(signal))
#' tss
#' stopifnot(tss == signal[1])
#'
#' # Neighboring "bonus" counts break such ties:
#'
#' # Look 1 basepair ahead to choose the maximum.
#' signal_ahead <- c(signal, GRanges("chr1:201:+", score = 1L))
#' signal_ahead
#' tss_ahead <- tss(signal_ahead, range(signal_ahead))
#' tss_ahead
#' stopifnot(tss_ahead == signal_ahead[2])
#'
#' # Look 1 basepair behind to choose the maximum.
#' signal_behind <- c(GRanges("chr1:99:+", score = 2L), signal_ahead)
#' signal_behind
#' tss_behind <- tss(signal_behind, range(signal_behind))
#' tss_behind
#' stopifnot(tss_behind == signal_behind[2])
#'
#' # Both, look ahead and look behind tie breaking.
#' signal_both <- sort(c(GRanges("chr1:199:+", score = 2L), signal_behind))
#' signal_both
#' tss_both <- tss(signal_both, range(signal_both))
#' tss_both
#' stopifnot(tss_both == signal_both[4])
#'
#' # Return Pairs with found region.
#' pairs <- tss(signal_both, range(signal_both), pairs = TRUE)
#' pairs
#' stopifnot(first(pairs) == tss_both)
#' stopifnot(second(pairs) == range(signal_both))
#' @export
tss <- function(signal, ranges, pairs = FALSE) {
    ## Find the highest peak within each region using 2 bp windows.  Tiling 2
    ## bp windows is inefficient, so instead use a rolling max in the C layer.
    stopifnot(is(ranges, "GRanges"))
    stopifnot(is(signal, "GRanges"))
    if (any(strand(signal) == "*"))
        stop("signal must be stranded")
    stopifnot(isSorted(signal))
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
    ## For the positive strand, break ties by preferring the first peak, and
    ## for the negative strand, break ties by preferring the last peak.  This
    ## choice is based on the ground-truth biological CA-motif at initiation
    ## sites; originally the negative strand also chose the first peak, but
    ## this caused the CA-motif (rather the reverse complement TG-motif) to not
    ## be nearly as enriched as seen on the positive strand, indicating that an
    ## identical signal peak within the same promoter on the negative strand
    ## was stealing precedence.
    prefer_last <- as.integer(strand(ranges[unique(queryHits(ol))]) == "-")
    .Call(
        C_tss, PACKAGE = "tsshmm", indices_peak, groups, indices_signal,
        starts_signal, scores_signal, prefer_last)
    if (! pairs) return(signal[indices_peak])
    indices_groups <- which(indices_signal %in% indices_peak)
    Pairs(signal[indices_peak], ranges[groups[indices_groups]])
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
