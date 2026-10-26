#' Compute the unweighted Minkowski aggregation of powered distances
#'
#' @param x Finite, non-negative matrix of powered distances.
#' @param degrees Minkowski order used to power the components.
#' @return One-column matrix containing the row-wise Minkowski norm.
#' @export
indexFrechet <- function(x, degrees = 2) {
  X <- as.matrix(x)
  if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L ||
      any(!is.finite(X)) || any(X < 0)) {
    stop("'x' must be a finite, non-negative numeric matrix.", call. = FALSE)
  }
  if (length(degrees) != 1L || !is.finite(degrees) || degrees < 1) {
    stop("'degrees' must be a finite number greater than or equal to 1.", call. = FALSE)
  }
  z <- matrix(rowSums(X)^(1 / degrees), ncol = 1L,
              dimnames = list(rownames(X), "Frechet.Index"))
  z
}
