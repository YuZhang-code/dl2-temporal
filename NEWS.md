# dl2 0.3.2

* Sensitivity analysis now compares Minkowski orders using the iterative
  algorithm in every scenario and reports iteration counts and convergence
  status.

# dl2 0.3.1

* Added `sensitivityDL2()` to compare rankings and top/bottom group membership
  across Minkowski orders and weighting methods.

# dl2 0.3.0

* Clarified that polarity is defined relative to the target construct, and added
  a vulnerability example in which higher unemployment increases the index.
* Restored the iterative DL2 procedure: initialize with an unweighted Frechet
  norm, order indicators by their correlation with the previous composite,
  update sequential correction factors, and stop using an Lp convergence
  distance with the same order as the composite norm.
* Retained permutation-averaged Peña–Trapero weighting as `weight_method = "pena"`.

# dl2 0.2.0

* Made polarity-aware normalization optional in `dl2()`.
* Changed Min–Max bounds to use the complete input sample, rather than computing
  separate annual extrema.
* Replaced MARS/earth variable importance and iterative smoothing with the
  permutation-averaged Peña–Trapero correction-factor procedure from the fuzzy
  implementation.
* Implemented DL2 aggregation as a weighted Minkowski norm and retained both
  raw correction factors and normalized aggregation weights in the result.
* Removed runtime dependencies on `earth`, `dplyr`, and `corrgram`.
