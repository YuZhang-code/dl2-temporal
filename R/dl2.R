#' Build a DL2 composite indicator
#'
#' The function orients indicators according to their polarity (optionally),
#' initializes the composite with an unweighted Minkowski distance, orders
#' indicators by their association with the previous composite, recalculates
#' sequential Peña–Trapero factors, and repeats until consecutive scores are
#' sufficiently close. The selected Minkowski order is used both for the score
#' and for the convergence distance. A separate permutation-averaged Peña
#' method is also available. For a time series or panel, pass all periods
#' together so normalization bounds and data-driven weights are common to all
#' periods within that sample.
#'
#' @param x Numeric data frame or matrix. Rows are observations and columns are
#'   indicators.
#' @param polarity Indicators whose higher raw values correspond to a higher
#'   value of the target construct, specified by column indices or names.
#'   For a vulnerability index, for example, unemployment belongs here when
#'   higher unemployment means greater vulnerability. Unlisted indicators are
#'   treated as negative-polarity: higher raw values correspond to a lower
#'   value of the target construct (for example, a protective indicator).
#'   If `NULL`, all indicators are treated as negative-polarity.
#' @param err Convergence tolerance for the iterative method, based on the
#'   normalized Lp distance between consecutive composite-score vectors.
#' @param iterations Maximum number of iterative updates.
#' @param degrees Minkowski order, at least 1. The default 2 gives the weighted
#'   Euclidean norm.
#' @param nfold,trace Deprecated MARS arguments; ignored.
#' @param use_normat If `TRUE`, apply polarity-aware Min–Max normalization
#'   using extrema across all rows in `x` (or supplied bounds). If `FALSE`, use
#'   the input values without normalization.
#' @param lambda_w_fixed Deprecated MARS smoothing argument; ignored.
#' @param custom_min,custom_max Optional fixed normalization bounds, one value
#'   per indicator. Supply both to apply previously estimated bounds.
#' @param reference Reference profile for the distance. Defaults to zero in the
#'   oriented scale; may be a named or positional numeric vector.
#' @param weights Optional fixed non-negative correction factors. When
#'   supplied, they are used directly and the data-driven iterative weighting
#'   step is bypassed. Otherwise, `dl2()` uses the iterative method.
#' @param weight_method The automatically estimated DL2 index uses the
#'   iterative sequential Peña–Trapero method. Use `iterations = 1` for one
#'   update. Permutation-averaged, equal, and subjective factors are available
#'   separately through [weightFactors()].
#' @param subjective_weights Retained for call compatibility; use
#'   [weightFactors()] for standalone subjective factor estimation.
#' @param n_permutations Retained for call compatibility; use
#'   [weightFactors()] for standalone permutation-based factor estimation.
#' @param epsilon Clipping tolerance for sequential R-squared estimates.
#' @param seed Reproducibility seed; `NULL` uses the current RNG stream.
#' @param corr_threshold,r2_warning_threshold Diagnostic thresholds.
#' @param verbose_warning Show weighting diagnostics.
#'
#' @return A list with the oriented data, fixed normalization bounds, raw
#'   correction factors, normalized aggregation weights, mean sequential R2,
#'   composite indicator, weighted components, reference profile, and
#'   diagnostics, iteration history, and weights at each iteration. The score
#'   has the direction of distance from `reference`.
#' @importFrom stats setNames
#' @importFrom utils tail
#' @export
dl2 <- function(x, polarity = NULL, err = 1e-4, iterations = 20L,
                degrees = 2, nfold = NULL, trace = NULL,
                use_normat = TRUE, lambda_w_fixed = NULL,
                custom_min = NULL, custom_max = NULL, reference = NULL,
                weights = NULL,
                weight_method = "iterative",
                subjective_weights = NULL, n_permutations = 200L,
                epsilon = 1e-6, seed = 123,
                corr_threshold = 0.90, r2_warning_threshold = 0.95,
                verbose_warning = TRUE) {
  if (!is.null(nfold) || !is.null(trace) || !is.null(lambda_w_fixed)) {
    warning("'nfold', 'trace', and 'lambda_w_fixed' are obsolete and ignored; MARS smoothing is no longer used.", call. = FALSE)
  }
  if (length(err) != 1L || !is.finite(err) || err < 0) {
    stop("'err' must be one finite non-negative convergence tolerance.", call. = FALSE)
  }
  if (length(iterations) != 1L || !is.finite(iterations) || iterations < 1 ||
      iterations != as.integer(iterations)) {
    stop("'iterations' must be a positive integer.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("'epsilon' must lie strictly between 0 and 1.", call. = FALSE)
  }

  X <- .dl2_matrix(x)
  raw_data <- as.data.frame(X, check.names = FALSE)
  if (use_normat) {
    normat <- normalization(raw_data, polarity = polarity,
                            custom_max = custom_max, custom_min = custom_min)
  } else {
    if (!is.null(custom_min) || !is.null(custom_max)) {
      stop("'custom_min' and 'custom_max' apply only when 'use_normat = TRUE'.", call. = FALSE)
    }
    if (!is.null(polarity)) warning("'polarity' is ignored because 'use_normat = FALSE'.", call. = FALSE)
    normat <- raw_data
  }
  Z <- .dl2_matrix(normat)
  vars <- colnames(Z)

  if (is.null(reference)) {
    reference <- setNames(rep(0, ncol(Z)), vars)
  } else {
    if (!is.numeric(reference) || length(reference) != ncol(Z) ||
        any(!is.finite(reference))) {
      stop("'reference' must be a finite numeric vector with one value per indicator.", call. = FALSE)
    }
    if (!is.null(names(reference))) {
      if (!setequal(names(reference), vars)) stop("Names of 'reference' must match indicator names.", call. = FALSE)
      reference <- reference[vars]
    } else {
      names(reference) <- vars
    }
  }

  score <- NULL
  wfit <- NULL
  variables_sort <- vars
  convergence_history <- numeric(0)
  composite_history <- NULL
  weights_by_iteration <- NULL

  if (is.null(weights)) {
    if (length(weight_method) != 1L || !identical(weight_method, "iterative")) {
      stop("'dl2()' uses the iterative weighting method. Set 'iterations = 1' for one update; use 'weightFactors()' for standalone alternative factors.", call. = FALSE)
    }
    weight_method <- "iterative"
    if (weight_method == "iterative") {
      # Keep the whole-sample coordinate scaling from normalization.R. The
      # iteration reorders indicators and updates their sequential factors.
      partial <- abs(sweep(Z, 2L, reference, "-"))
      dist.in <- partial^degrees
      previous <- as.numeric(indexFrechet(dist.in, degrees = degrees))
      composite_history <- matrix(previous, ncol = 1L,
                                  dimnames = list(NULL, "Frechet.initial"))
      weights_by_iteration <- matrix(numeric(0), nrow = 0L, ncol = ncol(Z),
                                     dimnames = list(NULL, vars))

      for (k in seq_len(as.integer(iterations))) {
        assoc <- suppressWarnings(vapply(seq_len(ncol(Z)), function(j) {
          stats::cor(partial[, j], previous)
        }, numeric(1)))
        assoc[!is.finite(assoc)] <- 0
        ord <- order(-abs(assoc), seq_along(assoc))
        variables_sort <- vars[ord]

        q <- setNames(numeric(ncol(Z)), vars)
        r2_seq <- setNames(numeric(ncol(Z)), vars)
        for (j in seq_along(ord)) {
          col_j <- ord[j]
          if (length(unique(partial[, col_j])) < 2L) {
            q[col_j] <- 0
            r2_seq[col_j] <- 1
          } else if (j == 1L) {
            q[col_j] <- 1
            r2_seq[col_j] <- 0
          } else {
            dat <- as.data.frame(partial[, c(ord[seq_len(j - 1L)], col_j), drop = FALSE])
            names(dat) <- make.unique(c(paste0("z", seq_len(j - 1L)), "y"))
            lm_fit <- tryCatch(stats::lm(y ~ ., data = dat), error = function(e) NULL)
            r2 <- if (is.null(lm_fit)) 0 else summary(lm_fit)$r.squared
            if (!is.finite(r2)) r2 <- 0
            r2 <- min(max(r2, 0), 1 - epsilon)
            r2_seq[col_j] <- r2
            q[col_j] <- 1 - r2
          }
        }
        if (!is.finite(sum(q)) || sum(q) <= 0) {
          stop("Iterative Pena-Trapero factors sum to zero; remove constant indicators or inspect the input.", call. = FALSE)
        }

        new_score <- as.numeric(calculateDP2(
          dist.in, q, degrees = degrees, iteration = k
        )$DP1)
        # Normalized Lp distance between consecutive score vectors. The same
        # p used in the composite is used in the convergence test.
        delta <- mean(abs(new_score - previous)^degrees)^(1 / degrees)
        convergence_history[k] <- delta
        weights_by_iteration <- rbind(weights_by_iteration, q / sum(q))
        rownames(weights_by_iteration)[k] <- paste0("iteration_", k)
        composite_history <- cbind(composite_history, new_score)
        colnames(composite_history)[k + 1L] <- paste0("DL2.", k)
        raw_factors <- q
        effective_weights <- q / sum(q)
        previous <- new_score
        if (delta <= err) break
      }

      score <- calculateDP2(dist.in, raw_factors, degrees = degrees,
                            iteration = length(convergence_history))
      cor_mat <- suppressWarnings(stats::cor(Z))
      high_pairs <- NULL
      if (ncol(Z) > 1L) {
        hi <- which(upper.tri(cor_mat) & is.finite(cor_mat) &
                      abs(cor_mat) >= corr_threshold, arr.ind = TRUE)
        if (nrow(hi)) high_pairs <- data.frame(
          var1 = vars[hi[, 1L]], var2 = vars[hi[, 2L]],
          correlation = cor_mat[hi], row.names = NULL
        )
      }
      high_r2 <- names(r2_seq)[is.finite(r2_seq) & r2_seq >= r2_warning_threshold]
      if (verbose_warning && length(high_pairs)) {
        warning("Highly correlated indicators detected; inspect 'high_correlation_pairs'.", call. = FALSE)
      }
      if (verbose_warning && length(high_r2)) {
        warning("High sequential R-squared values detected for: ",
                paste(high_r2, collapse = ", "), call. = FALSE)
      }
      wfit <- list(
        weights = raw_factors, correction_factors = raw_factors,
        aggregation_weights = effective_weights, R2 = r2_seq,
        raw_R2_by_permutation = NULL, method = "iterative",
        correlation_matrix = cor_mat, high_correlation_pairs = high_pairs,
        high_r2_variables = high_r2
      )
    } else {
      wfit <- weightFactors(
        Z, method = weight_method, subjective_weights = subjective_weights,
        n_permutations = n_permutations, epsilon = epsilon, seed = seed,
        corr_threshold = corr_threshold,
        r2_warning_threshold = r2_warning_threshold,
        verbose_warning = verbose_warning
      )
      raw_factors <- wfit$weights
      effective_weights <- wfit$aggregation_weights
      dist.in <- calculateDistance(Z, degrees = degrees, reference = reference)
      score <- calculateDP2(dist.in, raw_factors, degrees = degrees, iteration = 1L)
      composite_history <- as.matrix(score$DP1)
      colnames(composite_history) <- "DL2.1"
      weights_by_iteration <- matrix(effective_weights, nrow = 1L,
                                     dimnames = list("iteration_1", vars))
      convergence_history <- NA_real_
    }
  } else {
    if (!is.numeric(weights) || length(weights) != ncol(Z) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
      stop("'weights' must be a finite, non-negative vector with positive sum and one value per indicator.", call. = FALSE)
    }
    if (!is.null(names(weights))) {
      if (!setequal(names(weights), vars)) stop("Names of 'weights' must match indicator names.", call. = FALSE)
      weights <- weights[vars]
    }
    raw_factors <- setNames(as.numeric(weights), vars)
    effective_weights <- raw_factors / sum(raw_factors)
    wfit <- list(
      weights = raw_factors, correction_factors = raw_factors,
      aggregation_weights = effective_weights,
      R2 = NULL, correlation_matrix = suppressWarnings(stats::cor(Z)),
      high_correlation_pairs = NULL, high_r2_variables = character(),
      raw_R2_by_permutation = NULL, method = "user"
    )
    dist.in <- calculateDistance(Z, degrees = degrees, reference = reference)
    score <- calculateDP2(dist.in, raw_factors, degrees = degrees, iteration = 1L)
    composite_history <- as.matrix(score$DP1)
    colnames(composite_history) <- "DL2.1"
    weights_by_iteration <- matrix(effective_weights, nrow = 1L,
                                   dimnames = list("iteration_1", vars))
    convergence_history <- NA_real_
  }

  frechet <- indexFrechet(dist.in, degrees = degrees)
  bounds <- attr(normat, "bounds")
  invisible(list(
    normat = normat,
    use_normat = use_normat,
    normalization_bounds = bounds,
    Weights = effective_weights,
    aggregation_weights = effective_weights,
    correction_factors = raw_factors,
    R2 = wfit$R2,
    raw_R2_by_permutation = wfit$raw_R2_by_permutation,
    weight_method = wfit$method,
    derror = convergence_history,
    convergence_distance = if (length(convergence_history) &&
                               is.finite(tail(convergence_history, 1L)))
      tail(convergence_history, 1L) else NA_real_,
    composite_history = composite_history,
    weights_by_iteration = weights_by_iteration,
    frechet = frechet,
    corrgr = NULL,
    correlation_matrix = wfit$correlation_matrix,
    high_correlation_pairs = wfit$high_correlation_pairs,
    high_r2_variables = wfit$high_r2_variables,
    Comp.Indicator = score$DP1,
    Comp.Indicator_1 = score$DP1,
    iteration = if (wfit$method == "iterative") length(convergence_history) else 1L,
    variables_sort = variables_sort,
    dist.in = dist.in,
    dist.fin = score$mDP,
    reference = reference,
    base_delta = NA_real_
  ))
}
