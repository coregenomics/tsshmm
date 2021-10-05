context("encode")

test_that("encode checks inputs", {
    windows <- tile(ranges, width = 10)
    expect_silent(encode(signal, bg, windows))
    expect_error(encode(NULL, bg, windows), "signal")
    expect_error(encode(signal, NULL, windows), "bg")
    expect_error(encode(signal, bg, NULL), "windows")
})
