## Convert GRO-cap and GRO-seq files into GRangesList data object.

suppressPackageStartupMessages({
    library(BiocFileCache)
    library(GEOquery)
    library(rtracklayer)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(tibble)
    library(tidyr)
    library(usethis)
})

files_remote <- function(gse) {
    gds <- getGEO(gse)[[1]]

    pData(phenoData(gds)) %>%
        select(starts_with("supp")) %>%
        rownames_to_column("acc") %>%
        as_tibble(rownames) %>%
        mutate_if(is.factor, as.character) %>%
        gather(key = "type", value = "file", -acc) %>%
        select(-type) %>%
        filter(str_detect(file, "K562")) %>%
        pull(file)
}
files_remote <- unlist(map(c("GSE60453", "GSE60454"), files_remote))

cached_files <- function(files_remote, cache = getBFCOption("CACHE")) {
    bfc <- BiocFileCache(ask = FALSE, cache = cache)

    ## Assume the remote files do not change.  Simply check if they are present
    ## in the cache, otherwise download.
    tbl <- map_df(files_remote, bfcquery, x = bfc)
    files_to_add <- setdiff(files_remote, tbl$rname)
    ## Vectorized download.
    walk(files_to_add, bfcadd, x = bfc, rtype = "web")

    ## Get local paths.
    tbl <- map_df(files_remote, bfcquery, x = bfc)
    files_local <- bfcpath(bfc, tbl$rid)

    names(files_local) <- names(files_remote)
    files_local
}
files_local <- cached_files(files_remote)

core2014full <- as(map(files_local, import), "GRangesList")

## Set strand and treatment name.
assay <- str_extract(files_remote, "(GRO.{3})")
assay_type <- str_extract(files_remote, "([^_]+TAP)")
strand <- ifelse(str_detect(str_extract(files_remote, "(minus|plus)"),
                            "minus"), "-", "+")
trt <- str_c(assay, "_", assay_type)
trt[is.na(trt)] <- assay[is.na(assay_type)]
names(core2014full) <- trt
core2014full <- mendoapply(FUN = `strand<-`, core2014full, value = strand)
core2014full <- mendoapply(FUN = `score<-`, core2014full,
                           value = lapply(core2014full,
                                          function(x) abs(score(x))))

## Regroup and sort.
core2014full <- regroup(core2014full, trt)
core2014full <- sort(core2014full)

## > pryr::object_size(core2014full)
## 295,371,720 B

## Subset to bring down 64 MB rda file to about 200 kb.
gr <- range(unlist(core2014full, use.names = FALSE))
gr <- dropSeqlevels(gr, "chrM", "coarse")     # chrM is problematic.
target <- 0.2 / 64                            # 0.003125
actual <- min(lengths(gr)) / sum(lengths(gr)) # 0.005848092
subset <- gr[lengths(gr) == min(lengths(gr))]
core2014 <- endoapply(core2014full, subsetByOverlaps, subset)

## > pryr::object_size(core2014)
## 3,412,312 B

use_data(core2014, overwrite = TRUE)

## Don't save full data under data/ because it slows down load time!
save(core2014full,
     file = proj_path("data-raw", "core2014full.rda"),
     compress = formals(use_data)$compress,
     version = formals(use_data)$version)
