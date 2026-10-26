#' Assess sensitivity of a DL2 composite indicator
#'
#' Recalculates the index over combinations of Minkowski orders using the
#' iterative DL2 algorithm in every scenario, then compares each result with a
#' reference specification. Rank agreement and overlap among the most and
#' least vulnerable observations are
#' emphasized because score levels can change mechanically when the norm or
#' weights change.
#'
#' @param x Numeric matrix or data frame of indicators.
#' @param p_values Minkowski orders to compare. Every value must be at least 1.
#' @param weight_methods Must be `"iterative"`. Every scenario recalculates
#'   the composite until convergence or the maximum number of iterations.
#' @param reference_p Minkowski order defining the reference specification.
#' @param reference_method Weighting method defining the reference
#'   specification.
#' @param top_fraction Fraction of observations used to compare the top and
#'   bottom of the index ranking.
#' @param polarity,use_normat,custom_min,custom_max,reference,iterations,err,
#'   epsilon,seed,verbose_warning Arguments forwarded to [dl2()].
#' @param n_permutations Retained for compatibility. It has no effect because
#'   sensitivity scenarios use the iterative method, not permutation averaging.
#'
#' @return A list with a scenario comparison table, the reference index,
#'   indices from every scenario, and settings. The table reports Pearson score
#'   correlation, Spearman and Kendall rank correlations, normalized mean
#'   absolute score change, and top/bottom group overlap. High rank agreement
#'   and overlap indicate robustness to the tested choices; they do not
#'   establish construct or predictive validity. The table also reports the
#'   actual iteration count, convergence status, and final convergence
#'   distance for each scenario.
#' @export
sensitivityDL2 <- function(x, p_values = c(1, 2, 3),
                           weight_methods = "iterative",
                           reference_p = 2, reference_method = "iterative",
                           top_fraction = 0.20, polarity = NULL,
                           use_normat = TRUE, custom_min = NULL,
                           custom_max = NULL, reference = NULL,
                           iterations = 20L, err = 1e-4,
                           n_permutations = 200L, epsilon = 1e-6,
                           seed = 123, verbose_warning = FALSE) {
  X <- .dl2_matrix(x)
  allowed <- "iterative"
  if (!is.numeric(p_values) || !length(p_values) ||
      any(!is.finite(p_values)) || any(p_values < 1)) {
    stop("'p_values' must contain finite Minkowski orders greater than or equal to 1.", call. = FALSE)
  }
  if (!is.numeric(reference_p) || length(reference_p) != 1L ||
      !is.finite(reference_p) || reference_p < 1) {
    stop("'reference_p' must be one finite value greater than or equal to 1.", call. = FALSE)
  }
  if (!length(weight_methods) || any(!weight_methods %in% allowed) ||
      length(reference_method) != 1L || !reference_method %in% allowed) {
    stop("Sensitivity analysis requires the iterative method. Set 'iterations = 1' for one update.", call. = FALSE)
  }
  if (length(top_fraction) != 1L || !is.finite(top_fraction) ||
      top_fraction <= 0 || top_fraction > 0.5) {
    stop("'top_fraction' must lie in (0, 0.5].", call. = FALSE)
  }

  p_values <- unique(c(as.numeric(p_values), as.numeric(reference_p)))
  weight_methods <- unique(c(weight_methods, reference_method))
  scenarios <- expand.grid(
    p = p_values, weight_method = weight_methods,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  scenarios <- scenarios[order(scenarios$p, scenarios$weight_method), , drop = FALSE]
  fit_scenario <- function(p, method) {
    dl2(X, polarity = polarity, use_normat = use_normat,
        custom_min = custom_min, custom_max = custom_max,
        reference = reference, iterations = iterations, err = err,
        degrees = p, weight_method = method,
        n_permutations = n_permutations, epsilon = epsilon,
        seed = seed, verbose_warning = verbose_warning)
  }
  keys <- paste(scenarios$p, scenarios$weight_method, sep = "::")
  index_list <- setNames(vector("list", nrow(scenarios)), keys)
  fit_list <- setNames(vector("list", nrow(scenarios)), keys)
  for (i in seq_len(nrow(scenarios))) {
    fit <- fit_scenario(scenarios$p[i], scenarios$weight_method[i])
    index_list[[i]] <- as.numeric(fit$Comp.Indicator)
    fit_list[[i]] <- fit
  }
  ref_key <- paste(reference_p, reference_method, sep = "::")
  ref_score <- index_list[[ref_key]]
  n <- length(ref_score)
  n_group <- max(1L, ceiling(n * top_fraction))
  ref_order <- order(ref_score, seq_len(n))
  ref_bottom <- ref_order[seq_len(n_group)]
  ref_top <- tail(ref_order, n_group)
  ref_range <- diff(range(ref_score))

  metrics <- lapply(seq_len(nrow(scenarios)), function(i) {
    score <- index_list[[i]]
    rho <- suppressWarnings(stats::cor(ref_score, score, method = "spearman"))
    tau <- suppressWarnings(stats::cor(ref_score, score, method = "kendall"))
    pearson <- suppressWarnings(stats::cor(ref_score, score, method = "pearson"))
    scenario_order <- order(score, seq_len(n))
    bottom <- scenario_order[seq_len(n_group)]
    top <- tail(scenario_order, n_group)
    data.frame(
      p = scenarios$p[i], weight_method = scenarios$weight_method[i],
      is_reference = scenarios$p[i] == reference_p &&
        scenarios$weight_method[i] == reference_method,
      iterations_used = fit_list[[i]]$iteration,
      converged = is.finite(fit_list[[i]]$convergence_distance) &&
        fit_list[[i]]$convergence_distance <= err,
      final_convergence_distance = fit_list[[i]]$convergence_distance,
      pearson_score_correlation = pearson,
      spearman_rank_correlation = rho,
      kendall_rank_correlation = tau,
      normalized_mean_absolute_change = if (ref_range > 0)
        mean(abs(score - ref_score)) / ref_range else NA_real_,
      top_group_overlap = length(intersect(ref_top, top)) / n_group,
      bottom_group_overlap = length(intersect(ref_bottom, bottom)) / n_group,
      row.names = NULL
    )
  })
  names(index_list) <- keys
  list(
    comparison = do.call(rbind, metrics),
    reference_index = ref_score,
    indices = index_list,
    fits = fit_list,
    settings = list(reference_p = reference_p,
                    reference_method = reference_method,
                    top_fraction = top_fraction)
  )
}
