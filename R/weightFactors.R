#' Estimate Peña–Trapero correction factors
#'
#' For each random ordering of indicators, this function regresses each
#' indicator after the first on all preceding indicators. It averages each
#' indicator's sequential R-squared over permutations and returns the
#' redundancy correction factor `1 - mean(R2)`. These raw correction factors
#' are not constrained to sum to one; `dl2()` normalizes them only for the
#' weighted-norm aggregation.
#'
#' @param x Numeric matrix or data frame of indicators. Use the full reference
#'   sample (for a panel, all years) to obtain one set of weights for the
#'   entire series.
#' @param method Weighting method: `"pena"`, `"equal"`, or `"subjective"`.
#' @param subjective_weights Non-negative vector used when `method =
#'   "subjective"`.
#' @param n_permutations Positive number of random indicator orderings.
#' @param epsilon Positive clipping tolerance for sequential R-squared values.
#' @param seed Optional random seed. The caller's random-number state is
#'   restored after the calculation.
#' @param corr_threshold Threshold for reporting highly correlated pairs.
#' @param r2_warning_threshold Threshold for warning about high mean R-squared.
#' @param verbose_warning Show diagnostic warnings.
#'
#' @return A list containing raw `weights`/`correction_factors`, normalized
#'   `aggregation_weights`, mean `R2`, the correlation matrix, high-correlation
#'   pairs, high-R2 indicators, and the R2 values from each permutation.
#' @export
weightFactors <- function(x, method = c("pena", "equal", "subjective"),
                          subjective_weights = NULL, n_permutations = 200L,
                          epsilon = 1e-6, seed = 123,
                          corr_threshold = 0.90,
                          r2_warning_threshold = 0.95,
                          verbose_warning = TRUE) {
  method <- match.arg(method)
  X <- .dl2_matrix(x)
  n <- nrow(X)
  p <- ncol(X)
  vars <- colnames(X)

  if (length(n_permutations) != 1L || !is.finite(n_permutations) ||
      n_permutations < 1 || n_permutations != as.integer(n_permutations)) {
    stop("'n_permutations' must be a positive integer.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("'epsilon' must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (method == "equal") {
    raw <- setNames(rep(1, p), vars)
    R2 <- setNames(rep(0, p), vars)
    by_permutation <- matrix(0, nrow = 0, ncol = p,
                             dimnames = list(NULL, vars))
  } else if (method == "subjective") {
    if (is.null(subjective_weights) || length(subjective_weights) != p ||
        !is.numeric(subjective_weights) || any(!is.finite(subjective_weights)) ||
        any(subjective_weights < 0) || sum(subjective_weights) <= 0) {
      stop("'subjective_weights' must be a finite, non-negative vector with positive sum and one value per indicator.", call. = FALSE)
    }
    if (!is.null(names(subjective_weights))) {
      if (!setequal(names(subjective_weights), vars)) stop("Names of 'subjective_weights' must match indicator names.", call. = FALSE)
      subjective_weights <- subjective_weights[vars]
    }
    raw <- setNames(as.numeric(subjective_weights), vars)
    R2 <- setNames(rep(NA_real_, p), vars)
    by_permutation <- matrix(NA_real_, nrow = 0, ncol = p,
                             dimnames = list(NULL, vars))
  } else {
    if (n <= p) stop("Pena-Trapero weights require more observations than indicators to avoid saturated sequential regressions.", call. = FALSE)

    # Hold and restore the user's RNG state so a weight calculation does not
    # change randomness elsewhere in the user's analysis.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (!is.null(seed)) {
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }
    }, add = TRUE)
    if (!is.null(seed)) {
      if (length(seed) != 1L || !is.finite(seed)) stop("'seed' must be NULL or one finite number.", call. = FALSE)
      set.seed(as.integer(seed))
    }

    by_permutation <- matrix(0, nrow = as.integer(n_permutations), ncol = p,
                             dimnames = list(NULL, vars))
    constant <- vapply(seq_len(p), function(j) length(unique(X[, j])) < 2L,
                       logical(1))
    for (b in seq_len(as.integer(n_permutations))) {
      perm <- sample.int(p)
      ordered_r2 <- numeric(p)
      for (k in seq_len(p)) {
        j <- perm[k]
        if (constant[j]) {
          # A constant indicator has no discriminating information.
          ordered_r2[k] <- 1 - epsilon
        } else if (k == 1L) {
          ordered_r2[k] <- 0
        } else {
          dat <- as.data.frame(X[, c(perm[seq_len(k - 1L)], j), drop = FALSE])
          names(dat) <- make.unique(c(paste0("z", seq_len(k - 1L)), "y"))
          fit <- tryCatch(stats::lm(y ~ ., data = dat), error = function(e) NULL)
          r2 <- if (is.null(fit)) NA_real_ else summary(fit)$r.squared
          if (!is.finite(r2)) r2 <- 0
          ordered_r2[k] <- min(max(r2, 0), 1 - epsilon)
        }
      }
      by_permutation[b, perm] <- ordered_r2
    }

    R2 <- setNames(colMeans(by_permutation), vars)
    raw <- pmax(0, 1 - R2)
    names(raw) <- vars
    if (any(constant)) {
      raw[constant] <- 0
      R2[constant] <- 1
      by_permutation[, constant] <- 1
      if (verbose_warning) warning("Constant indicators receive zero correction factors.", call. = FALSE)
    }
    if (sum(raw) <= 0) {
      stop("All indicators have zero correction factors; no composite index can be formed.", call. = FALSE)
    }
  }

  corr_mat <- suppressWarnings(stats::cor(X))
  pairs <- NULL
  if (p >= 2L) {
    idx <- which(upper.tri(corr_mat) & is.finite(corr_mat) &
                   abs(corr_mat) >= corr_threshold, arr.ind = TRUE)
    if (nrow(idx) > 0L) {
      pairs <- data.frame(
        var1 = vars[idx[, 1L]], var2 = vars[idx[, 2L]],
        correlation = corr_mat[idx], row.names = NULL
      )
      pairs <- pairs[order(-abs(pairs$correlation)), , drop = FALSE]
      if (verbose_warning) warning("Highly correlated indicators detected; inspect 'high_correlation_pairs'.", call. = FALSE)
    }
  }

  high_r2 <- names(R2)[is.finite(R2) & R2 >= r2_warning_threshold]
  if (length(high_r2) && verbose_warning) {
    warning("High sequential R-squared values detected for: ",
            paste(high_r2, collapse = ", "), call. = FALSE)
  }

  w <- raw / sum(raw)
  names(w) <- vars
  list(
    weights = raw,
    aggregation_weights = w,
    correction_factors = raw,
    R2 = R2,
    correlation_matrix = corr_mat,
    high_correlation_pairs = pairs,
    high_r2_variables = high_r2,
    raw_R2_by_permutation = by_permutation,
    method = method
  )
}

#' @rdname weightFactors
#' @export
calculoFactoresPonderacion <- function(Z, ...) {
  weightFactors(Z, method = "pena", ...)
}
