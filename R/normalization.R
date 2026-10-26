#' Normalize indicators with common whole-sample bounds
#'
#' For each indicator, the extrema are computed over every row supplied to
#' `x`. In a time-series or panel application, pass all periods together so
#' that the same minimum and maximum are used in every year. Positive-polarity
#' indicators, whose higher raw values mean more of the target construct, are
#' oriented as `(x - min) / (max - min)`; negative-polarity
#' indicators are oriented as `(max - x) / (max - min)`.
#'
#' @param x Numeric matrix or data frame. Rows may represent units, years, or
#'   unit-year observations; columns are indicators.
#' @param polarity Indicators whose higher raw values correspond to more of
#'   the target construct, given as column indices or names. Every indicator
#'   not listed is treated as negative-polarity (higher raw values correspond
#'   to less of the target construct). If `NULL`, all indicators are treated
#'   as negative-polarity.
#' @param custom_max,custom_min Optional named or positional vectors of fixed
#'   bounds. Supply both to reuse bounds estimated on a reference sample.
#'
#' @return A data frame on the [0, 1] scale, with a `bounds` attribute
#'   containing the minima and maxima used for each column.
#' @export
normalization <- function(x, polarity = NULL, custom_max = NULL,
                          custom_min = NULL) {
  X <- .dl2_matrix(x)
  vars <- colnames(X)
  p <- ncol(X)

  if (xor(is.null(custom_min), is.null(custom_max))) {
    stop("Supply both 'custom_min' and 'custom_max', or neither.", call. = FALSE)
  }

  orient <- function(z, label) {
    if (is.null(z)) return(NULL)
    if (!is.numeric(z) || length(z) != p || any(!is.finite(z))) {
      stop(sprintf("'%s' must be a finite numeric vector with one value per indicator.", label), call. = FALSE)
    }
    if (!is.null(names(z))) {
      if (!setequal(names(z), vars)) stop(sprintf("Names of '%s' must match the indicator names.", label), call. = FALSE)
      z <- z[vars]
    }
    as.numeric(z)
  }

  lo <- orient(custom_min, "custom_min")
  hi <- orient(custom_max, "custom_max")
  if (is.null(lo)) {
    lo <- apply(X, 2L, min)
    hi <- apply(X, 2L, max)
  }
  if (any(hi < lo)) stop("Every maximum must be greater than or equal to its minimum.", call. = FALSE)

  if (is.null(polarity)) {
    positive <- integer(0)
  } else if (is.character(polarity)) {
    if (any(!polarity %in% vars)) stop("Unknown indicator name in 'polarity'.", call. = FALSE)
    positive <- match(polarity, vars)
  } else if (is.numeric(polarity) && all(is.finite(polarity)) &&
             all(polarity == as.integer(polarity)) &&
             all(polarity >= 1L & polarity <= p)) {
    positive <- as.integer(polarity)
  } else {
    stop("'polarity' must contain valid column indices or names.", call. = FALSE)
  }

  span <- hi - lo
  out <- matrix(0, nrow(X), p, dimnames = dimnames(X))
  variable <- span > 0
  for (j in which(variable)) {
    if (j %in% positive) {
      out[, j] <- (X[, j] - lo[j]) / span[j]
    } else {
      out[, j] <- (hi[j] - X[, j]) / span[j]
    }
  }
  # A constant indicator contains no variation and is mapped to zero.
  out <- as.data.frame(out, check.names = FALSE)
  attr(out, "bounds") <- data.frame(
    indicator = vars, minimum = as.numeric(lo), maximum = as.numeric(hi),
    polarity = ifelse(seq_len(p) %in% positive, "positive", "negative"),
    row.names = NULL
  )
  out
}

# Internal input validation shared by package functions.
.dl2_matrix <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  X <- as.matrix(x)
  if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L) {
    stop("'x' must have at least one row and one numeric indicator column.", call. = FALSE)
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  if (anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X))) {
    stop("Indicator names must be non-missing, non-empty, and unique.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("'x' must contain only finite values; resolve missing or infinite values before fitting DL2.", call. = FALSE)
  }
  X
}
