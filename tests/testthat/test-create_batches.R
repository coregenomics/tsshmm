context("create_batches helpers")

test_that("tile_with_rev is equivalent to endoapply(tile(...), rev)", {
    tiled <- tile(ranges, width = 10)
    expect_equal(tile_with_rev(ranges, 10, rev = FALSE), tiled)
    expect_equal(tile_with_rev(ranges, 10, rev = TRUE), endoapply(tiled, rev))
})
