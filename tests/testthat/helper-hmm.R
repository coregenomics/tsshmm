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
## Model trainined on full core2014 dataset with small B->P1 transition value
## that created an error with tsshmm:::prom_dist().  These parameters were
## generated from parameters(model_groseq) in 02-train.Rmd of:
## https://gitlab.com/coregenomics/tsshmmpaper
load(system.file("testdata", "parameters_model_groseq.rda",
                 package = "tsshmm", mustWork = TRUE))
