test_that("Mean of single position", {
  x <- rBigWig::fetch_region_means("testdata/test.bw", "1", 0, 1)
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
  expect_equal(x, 1)
})

test_that("Mean of two positions with same score", {
  x <- rBigWig::fetch_region_means("testdata/test.bw", "1", 1, 3)
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
  expect_equal(x, 1)
})

test_that("Mean of two positions at score-border", {
  x <- rBigWig::fetch_region_means("testdata/test.bw", "1", 99, 101)
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
  expect_equal(x, 1.5)
})

test_that("Mean at four positions at score-border with two bins", {
  x <- rBigWig::fetch_region_means("testdata/test.bw", "1", 99, 104, bins = 2)
  expect_is(x, "numeric")
  expect_equal(length(x), 2)
  expect_equal(x, c(1.5, 2))
})



test_that("Error if start == end", {
  expect_error( rBigWig::fetch_region_means("testdata/test.bw", "1", 1, 1) )
})

