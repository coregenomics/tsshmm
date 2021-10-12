context("model")

test_that("TSSHMM model initializes with bg_genebody toggle", {
    model_procap <- TSSHMM(bg_genebody = FALSE)
    expect_equivalent(dim(model_procap), c(7, 3))
    model_proseq <- TSSHMM(bg_genebody = TRUE)
    expect_equivalent(dim(model_proseq), c(8, 3))
})

test_that("TSSHMM model flags invalid bg_genebody input", {
    expect_error(TSSHMM(bg_genebody = 1.0))
    expect_error(TSSHMM(bg_genebody = -1))
    expect_error(TSSHMM(bg_genebody = 0))
    expect_error(TSSHMM(bg_genebody = "yes"))
    expect_error(TSSHMM(bg_genebody = "no"))
    expect_silent(TSSHMM(bg_genebody = TRUE))
    expect_silent(TSSHMM(bg_genebody = FALSE))
})

test_that("TSSHMM model is coerced to compact character strings", {
    model <- TSSHMM()
    lines <- as(model, "character")
    ## Without coercion you get string conversions like "c(0.9, 0.05, ...)"
    expect_false(any(grepl("c\\(", lines)))
    expect_false(any(grepl(",", lines)))
    ## No header; only numbers column joined on states dimension.
    dim <- dim(model)
    expect_equivalent(length(lines), dim["states"])
    ## Format should be "x.yz ... | x.yz ..."
    expect_true(all(grepl("\\|", lines)))
    chars <- table(strsplit(paste0(lines, collapse = " "), ""))
    symbols <- c(" ", ".", "|")
    numbers <- as.character(0:9)
    expect_equal(setdiff(names(chars),
                         c(symbols, numbers)),
                 vector("character"))
})

test_that("TSSHMM model supports both character coercion methods", {
    model <- TSSHMM()
    expect_equal(as(model, "character"), as.character(model))
})

test_that("TSSHMM show displays dimensions and parameters", {
    model <- TSSHMM()
    expect_output(show(model),
                  paste0("TSSHMM object with [0-9] hidden states ",
                         "and [0-9] emissions.*",
                         "Transition matrix:.*",
                         "Emission matrix:"))
})
