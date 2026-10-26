#' Aggregate powered distances with a weighted Minkowski norm
#'
#' @param x Non-negative powered-distance matrix, usually from
#'   [calculateDistance()].
#' @param xFactores Non-negative correction factors or weights, one per
#'   indicator. They are normalized to sum to one for aggregation.
#' @param degrees Same Minkowski order used to produce `x`.
#' @param iteration Integer suffix for the returned score column name.
#' @return A list with the composite score `DP1`, weighted powered components
#'   `mDP`, and normalized `weights`.
#' @export
calculateDP2 <- function(x, xFactores, degrees = 2, iteration = 1) {
  X <- as.matrix(x)
  if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L ||
      any(!is.finite(X)) || any(X < 0)) {
    stop("'x' must be a finite, non-negative numeric matrix.", call. = FALSE)
  }
  if (length(degrees) != 1L || !is.finite(degrees) || degrees < 1) {
    stop("'degrees' must be a finite number greater than or equal to 1.", call. = FALSE)
  }
  if (!is.numeric(xFactores) || length(xFactores) != ncol(X) ||
      any(!is.finite(xFactores)) || any(xFactores < 0) || sum(xFactores) <= 0) {
    stop("'xFactores' must be a finite, non-negative vector with positive sum and one value per column.", call. = FALSE)
  }
  w <- as.numeric(xFactores) / sum(xFactores)
  names(w) <- colnames(X)
  mDP <- sweep(X, 2L, w, "*")
  score <- rowSums(mDP)^(1 / degrees)
  DP1 <- matrix(score, ncol = 1L,
                dimnames = list(rownames(X), paste0("DL2.", iteration)))
  list(DP1 = DP1, mDP = mDP, weights = w)
}
