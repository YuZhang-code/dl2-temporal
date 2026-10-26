#' Compute powered distances from a reference profile
#'
#' @param x Numeric matrix or data frame of observations by indicators.
#' @param degrees Minkowski order, a finite number greater than or equal to 1.
#' @param reference Numeric reference vector, one value per indicator. Defaults
#'   to the zero vector.
#' @return A matrix containing `abs(x - reference)^degrees`.
#' @export
calculateDistance <- function(x, degrees = 2, reference = NULL) {
  X <- .dl2_matrix(x)
  if (length(degrees) != 1L || !is.finite(degrees) || degrees < 1) {
    stop("'degrees' must be a finite number greater than or equal to 1.", call. = FALSE)
  }
  if (is.null(reference)) reference <- rep(0, ncol(X))
  if (!is.numeric(reference) || length(reference) != ncol(X) ||
      any(!is.finite(reference))) {
    stop("'reference' must be a finite numeric vector with one value per indicator.", call. = FALSE)
  }
  if (!is.null(names(reference))) {
    if (!setequal(names(reference), colnames(X))) stop("Names of 'reference' must match indicator names.", call. = FALSE)
    reference <- reference[colnames(X)]
  }
  D <- abs(sweep(X, 2L, reference, "-"))^degrees
  colnames(D) <- colnames(X)
  D
}
