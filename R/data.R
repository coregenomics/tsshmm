#' Read counts of GRO-seq, and GRO-cap with and without TAP.
#'
#' Dataset of K562 cells from Core et al 2014.  The GRO-cap was differentially
#' treated with and without the critical 5'-cap removing enzyme, tobacco acid
#' phosphatase (TAP).
#'
#' @format GRangesList with 3 GRanges objects:
#' \describe{
#' \item{GROcap_wTAP:}{read counts, GRO-cap with TAP.}
#' \item{GROcap_noTAP:}{read counts, GRO-cap without TAP.}
#' \item{GROseq:}{read counts, GRO-seq.}
#' }
#' @source GSE60453 and GSE60454
#' @usage data(core2014)
"core2014"
