test_that("Pena factors are reproducible and preserve raw and normalized values", {
  x <- data.frame(a = c(1, 2, 4, 3, 8),
                  b = c(2, 1, 5, 2, 7),
                  c = c(8, 3, 6, 1, 5))
  w1 <- weightFactors(x, n_permutations = 40, seed = 19,
                      verbose_warning = FALSE)
  w2 <- weightFactors(x, n_permutations = 40, seed = 19,
                      verbose_warning = FALSE)
  expect_equal(w1$weights, 1 - w1$R2)
  expect_equal(w1$correction_factors, w1$weights)
  expect_equal(sum(w1$aggregation_weights), 1)
  expect_equal(w1$aggregation_weights, w2$aggregation_weights)
  expect_equal(w1$weights, w2$weights)
  expect_equal(w1$raw_R2_by_permutation, w2$raw_R2_by_permutation)
})

test_that("weighted Minkowski aggregation matches its definition", {
  powered <- matrix(c(9, 16), nrow = 1,
                    dimnames = list(NULL, c("a", "b")))
  ans <- calculateDP2(powered, c(a = 1, b = 3), degrees = 2)
  expect_equal(as.numeric(ans$DP1), sqrt(14.25))
  expect_equal(as.numeric(ans$weights), c(0.25, 0.75))
})

test_that("normalization can be disabled in the main algorithm", {
  x <- data.frame(a = c(2, 4, 5), b = c(4, 1, 3))
  ans <- suppressWarnings(dl2(x, use_normat = FALSE, weights = c(1, 1)))
  expect_equal(as.data.frame(ans$normat), x)
  expect_null(ans$normalization_bounds)
})

test_that("vulnerability polarity makes adverse movements raise the index", {
  x <- data.frame(unemployment = c(5, 10, 15), employment = c(80, 70, 60))
  ans <- dl2(x, polarity = "unemployment", weights = c(0.5, 0.5),
             use_normat = TRUE, verbose_warning = FALSE)
  expect_lt(ans$Comp.Indicator[1], ans$Comp.Indicator[2])
  expect_lt(ans$Comp.Indicator[2], ans$Comp.Indicator[3])
})

test_that("iterative convergence uses the selected Minkowski order", {
  x <- data.frame(a = c(1, 4, 2, 7, 5, 8),
                  b = c(8, 3, 6, 2, 7, 1),
                  c = c(2, 5, 3, 8, 4, 7))
  fit <- dl2(x, use_normat = FALSE, degrees = 1.5,
             iterations = 5, err = 0, verbose_warning = FALSE)
  hist <- fit$composite_history
  expect_equal(as.numeric(fit$Comp.Indicator), hist[, ncol(hist)])
  expect_length(fit$derror, fit$iteration)
  for (k in seq_along(fit$derror)) {
    expected <- mean(abs(hist[, k + 1L] - hist[, k])^1.5)^(1 / 1.5)
    expect_equal(fit$derror[k], expected)
  }
  expect_equal(nrow(fit$weights_by_iteration), fit$iteration)
})
