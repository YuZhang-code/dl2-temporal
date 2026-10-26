# A panel should be supplied as a complete sample, with one row per
# region-sector-year and one column per indicator. Do not normalize by year.
panel <- data.frame(
  unemployment_rate = c(8.0, 9.0, 12.0, 10.0),
  employment_rate = c(0.70, 0.68, 0.72, 0.74)
)

fit <- dl2(
  panel,
  # Vulnerability target: more unemployment raises vulnerability; more
  # employment lowers it, so employment is treated as negative-polarity.
  polarity = "unemployment_rate",
  use_normat = TRUE,
  degrees = 2,
  n_permutations = 500,
  seed = 123
)

fit$normalization_bounds
fit$correction_factors
fit$Weights
fit$Comp.Indicator
