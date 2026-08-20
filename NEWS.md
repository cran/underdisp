# underdisp 0.1.0

First public release. `underdisp` provides diagnostics and a unified family of
estimators for underdispersed count data (conditional variance below the
conditional mean), the case the Poisson and negative binomial defaults cannot
represent. The likelihood core is implemented in C++.

## Screening

* `ud_screen()`: at-risk dispersion diagnostic that benchmarks a fitted count
  against a zero-truncated Poisson under unit and period fixed effects, with a
  marginal statistic, a zero-truncated-Poisson at-risk statistic, and an
  over-conditioning guard. The at-risk threshold is available in two modes,
  `ztp_threshold = "calibrated"` (a fast n-dependent calibrated cutoff) and
  `ztp_threshold = "bootstrap"` (a parametric bootstrap of the fitted
  zero-truncated-Poisson null on the data's own design, `ztp_boot_B`), and an
  optional native COM-Poisson comparator (`run_comp`).

## Estimators

* `cpb()`: continuous parameter binomial (CPB) regression and its zero-truncated
  variant, with an interpretable observation-specific bound (the implied
  ceiling). C++ multistart estimation.
* `cpb_fe()`: high-dimensional fixed-effects CPB via a concentrated likelihood
  that matches the full-dummy likelihood to machine precision.
* `hurdle_cpb()` and `zi_cpb()`: two-part (participation logit + zero-truncated
  CPB intensity) and true zero-inflated mixture forms, each with fixed effects
  available in both equations.
* `gec()`: King's generalized event count / Katz-family regression, a single
  model whose dispersion is estimated freely and spans underdispersion,
  equidispersion, and overdispersion, with `gec_fe()`, `hurdle_gec()`, and
  `zi_gec()` variants.
* `count_reg()`, `hurdle_count()`, `zi_count()`: matched Poisson, negative
  binomial, and COM-Poisson families sharing one interface, with fixed effects
  (including two-way), offsets, and truncated/hurdle/zero-inflated forms, for
  like-for-like comparison against the CPB and GEC families.

## Inference

* Coefficient inference by a cold-multistart pairs (or cluster) bootstrap with
  percentile intervals; dispersion inference by a profile-likelihood interval
  for the shape parameter, which respects the feasibility boundary where the
  bootstrap distribution is skewed (`alpha_confint()`).
* Robust and clustered standard errors throughout, via a block bootstrap rather
  than a Hessian sandwich, since the CPB's parameter-dependent support makes the
  Hessian unreliable. Fixed-effects bootstraps use normal-approximation
  intervals, because replicates share the incidental-parameters bias and
  percentile intervals would be off-center.

## Bias-corrected fixed-effects estimation

* `bias_correct = "jackknife"` on the concentrated fixed-effects estimators
  (`cpb_fe()`, `gec_fe()`, and the `fe =` intensity of `hurdle_cpb()`): the
  split-panel jackknife (Dhaene & Jochmans 2015) removes the leading 1/T
  incidental-parameters bias of the dispersion parameter, validated by the
  package's fixed-effects bias Monte Carlo (alpha bias at T = 6/10/20/40 falls
  from -.147/-.096/-.056/-.033 to -.038/-.022/-.016/-.009). The correction
  carries a VALIDITY GATE for the method's time-homogeneity requirement: the
  time-split disagreement is compared against exchangeable unit-split placebo
  noise, and on panels where the two halves do not estimate a common parameter
  the correction is refused with an informative warning and the
  maximum-likelihood fit is returned. A refusal is diagnostic information
  about the panel, not an error.
* For separation in dummy-heavy participation equations, the hurdle
  factorizes: fit the participation stage with a dedicated bias-reduction
  package (`logistf`, `brglm2`) alongside the zero-truncated intensity from
  this package.

## Comparison and testing

* `compare_dispersion()` and `compare_models()`: information-criteria and proper
  scoring-rule comparison across Poisson, negative binomial, COM-Poisson, and
  CPB (and their hurdle and zero-inflated forms), with type-safe rejection of
  non-package or mismatched objects.
* `zi_test()`: boundary-corrected likelihood-ratio test for zero-inflation
  (the plain model nested in the mixture at the boundary), for both the CPB and
  GEC families.

## Quantities of interest

* `predict()`, `implied_ceiling()`, `alpha_confint()`, `irr()`, and
  `first_difference()`. Every `first_difference()` method in the package
  returns one unified `"ud_fd"` data frame (`component`, `from`, `to`, `diff`,
  `lower`, `upper`, `method`): single-equation fits report the mean, and the
  two-part/mixture models report the exact extensive/intensive decomposition
  (binary-stage probability, count-stage mean, marginal mean), with a `stage`
  argument on every two-part method to move the covariate in a single equation
  and a per-row `method` label recording the uncertainty source (stored
  bootstrap for `cpb`, refit bootstrap via `B=` for the CPB two-parts,
  delta-method where an analytic covariance exists). Unknown arguments error
  rather than being silently ignored.

## Distributions and simulators

* Density, distribution, quantile, and random-generation functions for the
  package's own families: `dcpb`/`pcpb`/`qcpb`/`rcpb`,
  `dgec`/`pgec`/`qgec`/`rgec`, `dcompois`/`pcompois`/`qcompois`/`rcompois`,
  plus the mixture simulators `rhurdle_cpb()` and `rzicpb()`. The standard
  zero-truncated, hurdle, and zero-inflated Poisson/negative-binomial forms are
  deliberately internal: CRAN already serves those families (e.g. `VGAM`), and
  the exported surface is reserved for what the package adds.
* `simulate()` methods for every model class (`cpb`, `cpb_fe`, `gec`, `gec_fe`,
  `count_reg`, `hurdle_count`, `zi_count`, `hurdle_cpb`, `hurdle_gec`, `zi_cpb`,
  `zi_gec`), with `fitted()` completed across the two-part and fixed-effects
  classes, so any fit plugs into `DHARMa`'s simulated-residual diagnostics via
  `DHARMa::createDHARMa()`.

## Panel tools and calibration

* `mundlak()`: correlated-random-effects device that augments a formula with
  unit means, model-agnostic across the estimators.
* `score()` and `cv_score()` (log score and ranked probability score, in sample
  and by cross-validation), `rootogram()` (Tukey hanging), and `pit_hist()`
  (non-randomized PIT histogram).

## Numerical robustness

The estimators are hardened for the difficult likelihoods this package targets
(near-boundary dispersion, zero-heavy and dummy-heavy panels):

* Multistart estimation with per-start restart-polishing, and, for the
  zero-inflated CPB, an expectation-maximization cross-check run alongside the
  direct maximum-likelihood multistart so the better optimum is kept, avoiding
  the degenerate zero-mass-collapse basin.
* Loud, estimator-backed fallbacks that degrade gracefully rather than returning
  `NA`: a Fisher-scoring fallback for the zero-truncated-Poisson screen when the
  default IRLS fit errors or silently diverges, and a ridge-regularized Fisher
  step in the bootstrap null on singular information matrices. Each fallback
  warns and reports the method actually used.
* Fitted-rate sanity guards that key on the pathology (diverged rates) rather
  than its correlates (iteration counts), and floored mixture likelihoods so a
  transient zero-probability iterate cannot produce a non-finite objective.
* Fit-rank–exact degrees of freedom in the at-risk screen and the mixture
  comparisons, so zero-heavy designs are neither over- nor under-penalized.

## Data and integration

* Bundled `peacekeeping` dataset (UN peacekeeping contributions, a
  state-year roster with structural zeros).
* `broom` (`tidy`/`glance`/`augment`) methods for all model classes, and
  `modelsummary` / `texreg` table support.
