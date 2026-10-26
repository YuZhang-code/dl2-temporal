\name{sensitivityDL2}
\alias{sensitivityDL2}
\title{Assess sensitivity of a DL2 composite indicator}
\usage{
sensitivityDL2(
  x,
  p_values = c(1, 2, 3),
  weight_methods = "iterative",
  reference_p = 2,
  reference_method = "iterative",
  top_fraction = 0.20,
  polarity = NULL,
  use_normat = TRUE,
  custom_min = NULL,
  custom_max = NULL,
  reference = NULL,
  iterations = 20L,
  err = 1e-04,
  n_permutations = 200L,
  epsilon = 1e-06,
  seed = 123,
  verbose_warning = FALSE
)
}
\arguments{
\item{x}{Numeric matrix or data frame of indicators.}

\item{p_values}{Minkowski orders to compare. Every value must be at least 1.}

\item{weight_methods}{Must be \code{"iterative"}. Every scenario recalculates
the composite until convergence or the maximum number of iterations.}

\item{reference_p}{Minkowski order defining the reference specification.}

\item{reference_method}{Weighting method defining the reference
specification.}

\item{top_fraction}{Fraction of observations used to compare the top and
bottom of the index ranking.}

\item{polarity, use_normat, custom_min, custom_max, reference, iterations,
err, epsilon, seed, verbose_warning}{Arguments forwarded to \code{dl2()}.}

\item{n_permutations}{Retained for compatibility. It has no effect because
sensitivity scenarios use the iterative method, not permutation averaging.}
}
\value{
A list with a scenario comparison table, the reference index, indices from
every scenario, and settings. The table reports Pearson score correlation,
Spearman and Kendall rank correlations, normalized mean absolute score change,
and top/bottom group overlap. High rank agreement and overlap indicate
robustness to the tested choices; they do not establish construct or predictive
validity. The table also reports the actual iteration count, convergence
status, and final convergence distance for each scenario.
}
\description{
Recalculates the index over combinations of Minkowski orders using the
iterative DL2 algorithm in every scenario, then compares each result with a
reference specification. Rank
agreement and overlap among the most and least vulnerable observations are
emphasized because score levels can change mechanically when the norm or
weights change.
}
\seealso{\code{dl2()}, \code{weightFactors()}}
\keyword{multivariate}
