test_that("normalization uses whole-sample extrema and respects polarity", {
  x <- data.frame(benefit = c(2, 4, 3), cost = c(10, 20, 15))
  z <- normalization(x, polarity = "benefit")
  expect_equal(z$benefit, c(0, 1, 0.5))
  expect_equal(z$cost, c(1, 0, 0.5))
  expect_equal(attr(z, "bounds")$minimum, c(2, 10))
  expect_equal(attr(z, "bounds")$maximum, c(4, 20))
})

test_that("polarity is oriented toward the target construct", {
  x <- data.frame(unemployment = c(5, 10, 15), employment = c(80, 70, 60))
  # The target is vulnerability: more unemployment and less employment mean
  # more vulnerability. Employment is therefore negative-polarity.
  z <- normalization(x, polarity = "unemployment")
  expect_equal(z$unemployment, c(0, 0.5, 1))
  expect_equal(z$employment, c(0, 0.5, 1))
  expect_gt(sum(z[3, ]), sum(z[1, ]))
})

test_that("fixed bounds can be reused and constant indicators are handled", {
  x <- data.frame(a = c(2, 4), b = c(7, 7))
  z <- normalization(x, polarity = "a", custom_min = c(a = 0, b = 7),
                     custom_max = c(a = 4, b = 7))
  expect_equal(z$a, c(0.5, 1))
  expect_equal(z$b, c(0, 0))
})
