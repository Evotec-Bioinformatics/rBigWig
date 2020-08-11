test_that("Retrieve single position", {
  x <- rBigWig::fetch_region("testdata/test.bw", "1", 0, 1)
  expect_is(x, "data.frame")
  expect_equal(ncol(x), 2)
  expect_true("Position" %in% colnames(x))
  expect_true("Score" %in% colnames(x))
  expect_equal(nrow(x), 1)
  expect_equal(x$Position, 0)
  expect_equal(x$Score, 1)
})


test_that("Retrieve two positions", {
  x <- rBigWig::fetch_region("testdata/test.bw", "1", 1, 3)
  expect_is(x, "data.frame")
  expect_equal(ncol(x), 2)
  expect_true("Position" %in% colnames(x))
  expect_true("Score" %in% colnames(x))
  expect_equal(nrow(x), 2)
  expect_equal(x$Position, c(1,2))
  expect_equal(x$Score, c(1,1))
})

test_that("Retrieve two positions at score-border", {
  x <- rBigWig::fetch_region("testdata/test.bw", "1", 99, 101)
  expect_is(x, "data.frame")
  expect_equal(ncol(x), 2)
  expect_true("Position" %in% colnames(x))
  expect_true("Score" %in% colnames(x))
  expect_equal(nrow(x), 2)
  expect_equal(x$Position, c(99,100))
  expect_equal(x$Score, c(1,2))
})

test_that("Error if start == end", {
  expect_error( rBigWig::fetch_region("testdata/test.bw", "1", 1, 1) )
})

