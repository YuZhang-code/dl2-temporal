\name{dl2-package}
\alias{dl2-package}
\alias{dl2}
\alias{normalization}
\alias{weightFactors}
\alias{calculoFactoresPonderacion}
\alias{calculateDistance}
\alias{calculateDP2}
\alias{indexFrechet}
\title{Iterative DL2 Composite Indicators with Global Temporal Normalization}
\description{
Constructs DL2 composite indicators from weighted Minkowski distances.
Normalization is optional. When used, extrema are calculated across every
row supplied to the function, so a panel or time series uses common bounds
across periods. The default weighting method iteratively orders indicators
by their absolute correlation with the preceding composite and recalculates
sequential Pena-Trapero factors. A separate permutation-averaged method is
available. Factors are normalized for weighted-norm aggregation.
}
\usage{
normalization(x, polarity = NULL, custom_max = NULL, custom_min = NULL)

weightFactors(x, method = c("pena", "equal", "subjective"),
  subjective_weights = NULL, n_permutations = 200L, epsilon = 1e-06,
  seed = 123, corr_threshold = 0.90, r2_warning_threshold = 0.95,
  verbose_warning = TRUE)

calculoFactoresPonderacion(Z, ...)

calculateDistance(x, degrees = 2, reference = NULL)

calculateDP2(x, xFactores, degrees = 2, iteration = 1)

indexFrechet(x, degrees = 2)

dl2(x, polarity = NULL, err = 1e-04, iterations = 20L, degrees = 2,
  nfold = NULL, trace = NULL, use_normat = TRUE, lambda_w_fixed = NULL,
  custom_min = NULL, custom_max = NULL, reference = NULL, weights = NULL,
  weight_method = "iterative",
  subjective_weights = NULL, n_permutations = 200L, epsilon = 1e-06,
  seed = 123, corr_threshold = 0.90, r2_warning_threshold = 0.95,
  verbose_warning = TRUE)
}
\arguments{
\item{x, Z}{Numeric matrix or data frame. Rows are observations and columns are indicators. For a panel, supply all periods together.}
\item{polarity}{Indicators whose higher raw values correspond to more of the target construct, specified by column indices or names. Unlisted indicators are negative-polarity: higher raw values correspond to less of the target construct. NULL treats all as negative. For a vulnerability index, unemployment is positive when higher unemployment means greater vulnerability; a protective indicator such as employment is negative.}
\item{custom_min, custom_max}{Optional fixed normalization bounds. Supply both to reuse bounds from a reference sample.}
\item{method}{Method used by code{weightFactors}: permutation-averaged Pena-Trapero, equal, or subjective.}
\item{weight_method}{The automatically estimated DL2 index uses the iterative sequential Pena-Trapero method. Use \code{iterations = 1} for one update. Permutation-averaged, equal, and subjective factors are available separately through \code{weightFactors()}.}
\item{subjective_weights}{Retained for call compatibility; use \code{weightFactors()} for standalone subjective factor estimation.}
\item{n_permutations}{Retained for call compatibility; use \code{weightFactors()} for standalone permutation-based factor estimation.}
\item{epsilon}{Clipping tolerance for sequential R-squared values.}
\item{seed}{Optional random seed; the user's RNG state is restored.}
\item{corr_threshold}{Absolute correlation threshold for diagnostics.}
\item{r2_warning_threshold}{Mean sequential R-squared threshold for diagnostics.}
\item{verbose_warning}{Whether to issue diagnostic warnings.}
\item{degrees}{Minkowski order, a finite number greater than or equal to 1.}
\item{reference}{Finite reference vector, one value per indicator; defaults to zero.}
\item{xFactores}{Finite non-negative correction factors or weights, one per indicator.}
\item{iteration}{Integer suffix for the resulting score column name.}
\item{err}{Convergence tolerance, calculated as the normalized Lp distance between consecutive composite-score vectors.}
\item{iterations}{Maximum number of iterative updates.}
\item{nfold, trace, lambda_w_fixed}{Deprecated MARS-era arguments, retained for call compatibility and ignored.}
\item{use_normat}{If TRUE, apply polarity-aware Min-Max normalization over all rows in x. If FALSE, use values as supplied.}
\item{weights}{Optional fixed non-negative correction factors. When supplied, they are used directly and the data-driven iterative weighting step is bypassed. Otherwise, \code{dl2()} uses the iterative method.}
\item{...}{Additional arguments passed by the Spanish alias to \code{weightFactors}.}
}
\details{
The iterative method initializes the score with the unweighted Frechet norm.
At each step, indicators are ordered by their absolute correlation with the
previous score. The first factor is one; each subsequent raw factor is one
minus the R-squared from regressing that indicator on all preceding
indicators in the current order. Factors are normalized for aggregation.
The score is recalculated until the normalized Lp distance between
consecutive score vectors is below err or the iteration limit is reached.
The separate permutation-averaged method averages sequential R-squared
values over random permutations and returns raw factors
\eqn{1 - \overline{R^2_j}}.

The score is \eqn{(\sum_j w_j |z_{ij}-b_j|^p)^{1/p}}. After optional
polarity-aware normalization, a larger value indicates greater distance from
the specified reference profile. Polarity is defined relative to the target
construct. For a vulnerability index, positive polarity means that higher
raw values indicate greater vulnerability; negative polarity means that
higher raw values indicate less vulnerability. When normalization is
disabled, the user is responsible for the scale, direction, and meaning of
the reference profile.
}
\value{
\code{normalization} returns a normalized data frame with the extrema in a
\code{bounds} attribute. \code{weightFactors} returns raw factors in
\code{weights} and \code{correction_factors}, normalized aggregation weights
in \code{aggregation_weights}, mean R-squared values, and diagnostics.
Distance helpers return matrices or lists of weighted components. \code{dl2} returns a list
with the data used, bounds, weights, correction factors, composite scores,
reference profile, and diagnostics.
}
\examples{
\dontrun{
x <- data.frame(unemployment = c(8, 10, 13), employment = c(70, 68, 64))
fit <- dl2(x, polarity = "unemployment", use_normat = TRUE, seed = 42)
fit$Comp.Indicator
}
}
\seealso{\code{\link{dl2}}, \code{\link{normalization}}, \code{\link{weightFactors}}}
