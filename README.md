# dl2 0.3.2

An R implementation of a DL2 composite indicator based on weighted Minkowski
distances, optional whole-sample normalization, and iterative sequential
Peña–Trapero redundancy correction factors.

## What changed from 0.1.0

- Min–Max normalization is optional (`use_normat = TRUE/FALSE`).
- If normalization is used, each indicator's minimum and maximum are computed
  across **all rows supplied**. For a panel, pass the full time series together
  to retain one scale across years. Fixed bounds can be reused with
  `custom_min` and `custom_max`.
- MARS/`earth` spline variable importance has been removed. The default
  `weight_method = "iterative"` follows the original iterative DL2 structure:
  initialize with an unweighted Frechet norm, order indicators by their
  absolute correlation with the previous composite, calculate sequential
  regression factors, and update the composite until its Lp change is below
  `err` or `iterations` is reached.
- `weightFactors(method = "pena")` retains the permutation-averaged procedure
  as a standalone factor estimator. `dl2()` itself uses the iterative
  algorithm; set `iterations = 1` to perform one update.
- Raw correction factors remain available in the output. They are normalized
  to sum to one only for the weighted-norm aggregation.
- The runtime package has no dependencies beyond base R.

## Formula

In the iterative method, at iteration $k$, indicators are ordered by their
absolute correlation with $I^{(k-1)}$. In that order, the first raw factor is
one and each subsequent factor is $q_j^{(k)}=1-R_{j,k}^2$, where $R_{j,k}^2$
comes from regressing the indicator on all preceding indicators. The raw
factors and normalized aggregation weights are

\[
q_j^{(k)} = 1 - R_{j,k}^2, \qquad
w_j^{(k)} = \frac{q_j^{(k)}}{\sum_l q_l^{(k)}}.
\]

For observation \(i\), the score is

\[
I_i^{(k)} = \left(\sum_j w_j^{(k)} |z_{ij} - b_j|^p\right)^{1/p},
\]

where \(z_{ij}\) is the input after optional polarity-aware normalization,
\(b\) is the reference profile, and \(p\geq 1\) is `degrees`. A positive
polarity means that a higher raw value represents more of the target construct;
a negative polarity means that it represents less. Thus, for a vulnerability
index, unemployment is positive when higher unemployment means greater
vulnerability, while employment or financial strength is negative when higher
values indicate greater resilience. With this orientation and the default zero
reference, larger scores represent greater distance in the direction of the
target construct (holding weights and other indicators fixed). Convergence is
assessed using the same order $p$:
\[
\Delta_k = \left(\frac{1}{n}\sum_i
  |I_i^{(k)}-I_i^{(k-1)}|^p\right)^{1/p}.
\]
The loop stops when \(\Delta_k\leq\texttt{err}\), or when it reaches
`iterations`.

## Example with time series

```r
panel <- data.frame(
  year = 2019:2022,
  unemployment = c(12.0, 13.5, 15.0, 14.0), # higher -> more vulnerability
  employment = c(0.62, 0.65, 0.60, 0.69)     # higher -> more resilience
)

# Keep the complete series together. The two extrema for each indicator are
# calculated once over all four observations, not separately by year.
fit <- dl2(
  panel[c("unemployment", "employment")],
  polarity = "unemployment",
  use_normat = TRUE,
  degrees = 2,
  iterations = 20,
  err = 1e-4
)

fit$Comp.Indicator
fit$correction_factors
fit$Weights
fit$normalization_bounds
```

Here the index measures vulnerability: unemployment has positive polarity,
whereas employment is unlisted and therefore has negative polarity. The
normalized unemployment component rises with unemployment; the normalized
employment component falls as employment rises. If instead the target were
resilience, these polarity assignments would be reversed.

To skip normalization, set `use_normat = FALSE`. In that case, check that
indicators have a meaningful common scale and orientation, and supply
`reference` explicitly when zero is not the intended baseline. If subsequent
data must be scored against a fixed historical reference, reuse the fitted
minimum and maximum values and the same weights rather than re-estimating them:

```r
bounds <- fit$normalization_bounds
new_fit <- dl2(
  new_panel[c("unemployment", "employment")],
  polarity = "unemployment",
  custom_min = setNames(bounds$minimum, bounds$indicator),
  custom_max = setNames(bounds$maximum, bounds$indicator),
  weights = fit$Weights,
  use_normat = TRUE
)
```

All input indicators must be finite; resolve missing and infinite values before
calling the package functions.

## Weight estimation alone

```r
weightFactors(
  panel[c("unemployment", "employment")],
  method = "pena", n_permutations = 500, seed = 123
)
```

The returned `weights` and `correction_factors` are the raw factors from the
selected standalone method; `aggregation_weights` contains their normalized
version. A constant indicator is assigned a zero correction factor because it
does not discriminate across observations. These standalone factors are not
used by `dl2()` unless supplied explicitly through its `weights` argument.

## Sensitivity and robustness

Convergence only shows that successive iterations are close. It does not show
that the index is insensitive to the selected norm or weighting rule. Compare
reasonable alternatives with `sensitivityDL2()`:

```r
robust <- sensitivityDL2(
  panel[c("unemployment", "employment")],
  polarity = "unemployment",
  p_values = c(1, 2, 3),
  reference_p = 2,
  reference_method = "iterative",
  iterations = 20,
  n_permutations = 500,
  seed = 123
)

robust$comparison
```

The comparison reports score correlations, Spearman and Kendall rank
correlations, mean absolute score changes relative to the reference index's
range, and the share of units retained in the highest and lowest `top_fraction`
of the ranking. Sensitivity scenarios all use the iterative DL2 method; the
table also reports how many updates each fit used and whether it reached the
convergence tolerance. Set `iterations = 1` to perform exactly one update.
Rank correlations and group overlap are the main robustness checks because
changing the norm can change score scale. This analysis assesses sensitivity
to the norm order; it does not establish construct validity, predictive
validity, or causality. Those require external benchmarks, outcome data, or a
separate identification strategy.
