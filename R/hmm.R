#' @importFrom BiocGenerics relist
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges findOverlaps tile
#' @importFrom S4Vectors endoapply mcols `mcols<-` to
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
hmm <- function(signal, bg, ranges) {
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
    windows <- tile(ranges, width = 10)
    observations <- encode(signal, bg, windows[has_signal])
    is_promoter <- endoapply(observations, viterbi)
    is_promoter <- is_promoter > 0

    reduce(unlist(windows[has_signal])[unlist(is_promoter)])
}

scoreOverlaps <- function(gr, reads) {
    ol <- findOverlaps(gr, reads)
    mcols(ol)$score <- score(reads[to(ol)])
    df <- as(ol, "DataFrame")
    df_min <- aggregate(score ~ queryHits, df, min)
    df_max <- aggregate(score ~ queryHits, df, max)
    score <- vector("integer", length(gr))
    score[df_min$queryHits] <- ifelse(
        abs(df_min$score) > abs(df_max$score),
        df_min$score, df_max$score)
    score
}

scoreOverlapWindows <- function(grl, reads)
    relist(scoreOverlaps(unlist(grl), reads), grl)
