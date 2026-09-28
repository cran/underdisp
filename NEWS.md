# underdisp 0.1.1

## New families

* `count_reg()`, `hurdle_count()`, and `zi_count()` gain four underdispersed
  families through the same interface, each inheriting fixed effects,
  zero-truncation, hurdle and zero-inflated forms, offsets, weights, the
  sandwiches, proper scores, DHARMa simulation, and the QoI methods:
  `"mpcmp"` (the COM-Poisson in Huang's 2017 mean parameterization, with the
  rate solving E(Y) = mu found in C++), `"genpois"` (the Consul-Jain generalized
  Poisson with constant dispersion lambda, renormalized on its finite support
  when lambda < 0), `"gammacount"` (Winkelmann's 1995 renewal-process model
  with gamma waiting times), and `"doublepois"` (Efron's 1986 double Poisson
  with the exact normalizing constant). Every family's pmf was checked against
  an independent implementation (rmutil, gamlss.dist, VGAM) to machine
  precision, and the generalized-Poisson and mean-parameterized COM-Poisson
  log-likelihoods reproduce glmmTMB's `genpois` and `compois` fits exactly.
* `dgammacount()`, `ddoublepois()`, `dgenpois()` and their p/q/r companions;
  the COM-Poisson functions accept `mu =` for the mean parameterization.

## New diagnostics

* `dispersion_test()`: the likelihood-ratio test of a fitted family's
  dispersion parameter against its Poisson value on the fitted design (fixed
  effects included, the unit effects profiled out of the Poisson null), with
  two-sided or directional alternatives by the signed root and the boundary
  mixture where the Poisson value is on the boundary (CPB, negative binomial);
  and the Cameron-Trivedi auxiliary-regression test for a Poisson fit. The
  p-value is a parametric bootstrap under the fitted Poisson (`B`, `cores`)
  whenever the Poisson value is on the boundary of the family's parameter space
  (CPB, negative binomial) or the fit has unit fixed effects; where it is
  interior, the asymptotic p-value is used unless the first-order bias from the
  null's estimated mean parameters, p / sqrt(2n) standard deviations, pushes a
  5% test's size above 6%. `summary()` of a CPB fit flags its asymptotic p-value.
* `ud_screen()`: in `ztp_threshold = "bootstrap"` mode the NB-vs-Poisson p-value
  is a parametric bootstrap under the fitted Poisson (`nb_boot`) when the mean
  model is within the comparator gates; the fast mode's asymptotic value is
  labeled conservative. The bootstrap restores the random-number state, so the
  at-risk threshold and every other result are unchanged.
* `dispersion_profile()`: the empirical conditional variance-to-mean ratio by
  bins of the fitted mean against the ratio each fitted family implies, with a
  `plot()` method, so the mechanism behind an underdispersed outcome (a hard
  ceiling, regular event timing, a soft tail) can be read off the data.
* `calibrate_alpha()`: an opt-in calibration of the interval for the CPB
  dispersion parameter. The support of the CPB moves with the parameters
  through the log-ceiling coefficients, which behave like endpoint parameters,
  and the estimate of `alpha` inherits a bias toward zero from them; the
  first-order (chi-square) profile interval therefore covers about 0.90 to
  0.94 in simulations from the CPB, with its misses on the upper side.
  `calibrate_alpha(fit, B, cores, seed)` stores the parametric-bootstrap
  distribution of the signed root of the profile likelihood ratio in the fit,
  and `alpha_confint()`, `confint()`, `implied_ceiling()` and `summary()` then
  report the interval cut at its quantiles (0.92 to 0.96 in the same
  simulations, misses balanced). The responses are drawn once on the calling
  process and the refits use no random numbers, so the result does not depend
  on `cores`. It is refused, with the reason, for a fit at the `max.support`
  guard, weights that are not whole numbers, and a fit saved by 0.1.0. With few
  distinct covariate patterns the first-order interval already covers at its
  level and the calibration is not needed. Without it the accessors report the
  first-order interval, and `attr(, "method")` and `summary()` say which one is
  shown.

## Weights, parallel bootstraps, and column scaling

* `weights =` (frequency weights) on every estimator; a weight of w is
  equivalent to w copies of the row in the likelihood, the information, the
  score products, and the bootstrap, and `nobs()` returns the weight total.
  In the parametric bootstraps (`dispersion_test()`, `calibrate_alpha()`) a row
  of weight w contributes w independent simulated responses; weights that are
  not whole numbers have no such replicate data, so `dispersion_test()` then
  reports the asymptotic p-value with a note.
* `cores =` on every bootstrap, the parametric-bootstrap screen threshold,
  `cv_score()`, and the refit first differences: replicates run on a socket
  cluster seeded from the calling session.
* Every optimizer works on covariates scaled to unit standard deviation and
  maps the coefficients back, and the Nelder-Mead fits restart from their own
  solution until the improvement stops, so a fit no longer depends on the
  units of a covariate. The numerical Hessians and score matrices of the count
  families are differenced in the same scaled parameterization.

## The CPB mean

* `fitted()`, `predict(type = "response")`, `residuals()`, and
  `first_difference(quantity = "mean")` for `cpb`, `cpb_fe`, and the hurdle-CPB
  intensity now report the exact mean of the fitted distribution (a finite sum
  of the renormalized pmf, and the conditional mean E(Y | Y >= 1) for
  zero-truncated fits), as `gec()` and the matched families already did. The
  rate parameter exp(x'b) equals the mean only when lambda/(1-alpha) is an
  integer; it remains available as `$rate` and `predict(type = "rate")`.
  The same holds for the generalized event count: `delta` is the Katz
  dispersion parameter, which is the variance-to-mean ratio on an unbounded
  support, and `exp(x'b)` is the rate of the recursion.
* Zero-truncated `count_reg()` fits report E(Y | Y >= 1).

## Methods

* `confint()` for every model class (bootstrap percentile, normal-approximation
  from the bootstrap, or Wald from the analytic covariance, one rule per
  inference type); `vcov()` for the hurdle and zero-inflated classes;
  `implied_ceiling()` for `cpb_fe`; `irr()` and `first_difference()` on a `cpb`
  fit without a bootstrap return point estimates with `method = "none"`;
  `vcov()` returns `NULL` rather than an error when no inference was requested.
* `irr()` labels the binary stage `"participation"` (hurdle) or `"inflation"`
  (zero-inflated), which point in opposite directions.
* texreg `extract()` methods for `count_reg`, `hurdle_count`, `zi_count`; the
  methods for every class now attach to texreg's generic (in 0.1.0 the
  registration silently failed and tables fell back to broom), and the
  count families label their dispersion row by the parameter.
* `tidy()` for the two-part classes prefixes each term by its component, so
  `modelsummary()` lays the table out without a `shape` argument.
* `predict()`, `score()`, and the QoI methods rebuild a `newdata` model matrix
  in the fit's own factor coding (stored levels and contrasts) and align it to
  the coefficients by name; a factor whose reference level was set through its
  contrasts attribute (as in AER's `NMES1988`) is no longer re-coded on
  `newdata`, and a coding mismatch is refused instead of multiplied by position.
* A formula `offset()` term and a pscl-style `y ~ x | z` formula are refused
  with a message naming `offset =`, `participation =`, and `zero =`; rows whose
  unit identifier is missing are dropped like other incomplete rows; the
  screen's generalized-Poisson comparator is the package's own family, so it
  can report underdispersion, and it is gated by `comp_max_par`/`comp_max_n`
  like the COM-Poisson (a full fit, not a diagnostic, on a dummy-heavy design); `print.gec()` labels the dispersion direction by
  the likelihood-ratio test; `count_reg(fe =)` hides the absorbed dummies and
  warns when `cluster` is the fixed-effects variable; standard errors from a
  Hessian that is not positive definite are `NA` with a warning; the
  generalized Poisson's support is capped for negative dispersion; `score()`
  honours frequency weights; the auxiliary test is weighted; the parallel
  driver draws its cluster seed eagerly, so consecutive `cores > 1` bootstraps
  differ.
* `converged` is stored by every class and reported by `print()`/`summary()`
  when it is `FALSE`; a matched-family fit whose optimizer never left the
  infeasible plateau (a fitted mean outside `[1e-12, 1e8]`) is reported as not
  converged with a warning, and its log-likelihood is `-Inf`, instead of
  returning the sentinel as a likelihood.
* The predictive distribution of a zero-truncated `gec()` or `count_reg()`
  fit, as `score()`, `rootogram()`, `pit_hist()`, and `cv_score()` use it, is
  the fit's own zero-truncated pmf normalized over the full support (in 0.1.0
  the Katz pmf was scored untruncated and the count families renormalized over
  the observed range), so the in-sample log score equals `-logLik/n` for every
  class; `fitted()` on a `hurdle_gec` is the exact conditional mean that
  `predict()` returns.
* `cv_score()`, `rootogram()`, `pit_hist()`, and `dispersion_profile()` apply
  frequency weights as `score()` does, so a row of weight w counts as w copies
  of itself: the rootogram's observed and expected frequencies, the PIT
  histogram, the dispersion profile's bins (equal weight, rows kept whole) and
  bin means, and, through a new `weights` argument, the mean held-out scores of
  `cv_score()` and `score(newdata =)`. `rootogram()` no longer fails when `kmax`
  exceeds the largest count.
* `pit_hist()` read the lower bound F(y - 1) out of row order wherever a zero
  count preceded a positive one, so its histogram was wrong for data with
  zeros; the bound is now looked up row by row.
* `ud_screen()` reports the NB-vs-Poisson likelihood ratio at its boundary
  value 0 when `glm.nb()` stops short of the Poisson limit, and prints "at the
  Poisson boundary" there rather than an overdispersion reading.
* The zero-truncated CPB at an implied ceiling below 1 is the point mass at 1
  in the pmf, the moments, the first differences, and the simulator alike (the
  zero-truncated distribution for every ceiling in [1, 2)), so the hurdle's
  predictions, scores, and rootograms are finite at rows whose covariates put
  the rate below a unit's feasibility floor.
* One rule names the coefficients of the two-part classes everywhere:
  `participation:`/`intensity:` for the hurdles and `count:`/`zero:` for the
  zero-inflated models, in `vcov()`, `confint()`, `tidy()`, and the texreg
  tables, so a term read off a table indexes the covariance and the interval.
* A vector `cluster` is checked for its length and for missing values by every
  estimator before it is used (the fixed-effects and two-part estimators
  recycled or padded a short vector in 0.1.0), and a row dropped for a missing
  unit identifier drops its entry of a vector `offset`, `weights`, or
  `cluster` with it; `hurdle_count()` selects the cluster-robust covariance
  when `cluster` is supplied, as its siblings do, and accepts a vector
  `cluster` aligned to the data; `zi_cpb()` and `zi_gec()` refuse a
  non-integer response instead of flooring it; `zi_cpb()`, `zi_count()`, and
  `zi_gec()` work on scaled covariates with a Nelder-Mead polish like every
  other estimator (the unscaled joint fit of `zi_count()` could stop at its
  starting values: on lme4's `grouseticks` it stopped 4.4 log-likelihood units
  below pscl's optimum); `summary.zi_count()` prints the
  stored inflation standard errors; every bootstrapping summary states how
  many replicates converged; `count_reg()` and `zi_count()` store
  `$fitted.values` on the mean scale and the natural parameter as `$mu`.
* The parallel workers of `cores > 1` attach the package, so a `cv_score()`
  fit function that names `cpb()` or `count_reg()` unqualified runs on them.
* The analytic, robust, and cluster covariances of the matched families are
  differenced at an absolute step of 1e-3 in the scaled parameterization
  (numDeriv's default relative step left the family's support at large fitted
  means and returned `NA` standard errors).
* `predict()` with `newdata` returns `NA` for a row whose covariate or offset is
  missing, as `predict.glm()` does, instead of failing inside the distribution
  routines. A column holding only missing values (which R reads as logical) is
  treated as a missing value of the fitted variable's type, and a variable
  supplied with a different type than in the fitting data is refused with a
  message that says so.

## Estimation and numerics

* The CPB profile-likelihood interval for `alpha` is traced by continuation
  with a feasibility-repaired start at each `alpha`, and the profile at each
  `alpha` is maximized with the fit's own multistart (BFGS, Nelder-Mead
  polish, perturbed restarts), so the lower limit is the re-maximized
  profile's crossing point; the interval is one-sided only when the profile
  has not dropped to the cut by `alpha = 0.005`. Each limit is the outermost
  crossing of its cut (the profile is a saw-tooth, so the set above a cut need
  not be connected), and a cut deeper than the stored trace, as at a level
  above 0.995, re-runs the trace instead of returning the parameter bounds.
* `rcpb(truncated = TRUE)` warns when a rate's ceiling `lambda / (1 - alpha)` is
  below 1: the zero-truncated CPB does not exist there, and the draw of 1 it
  returns has probability zero under the stated parameters.
* `max.support` defaults to `max(500, 10 * max(y))`; a fit whose largest
  implied ceiling reaches the guard warns and records `$support_binding`; a
  fit whose every start is infeasible errors instead of returning the starting
  values; a CPB fit that the guard stops short of its Poisson limit reports the
  likelihood-ratio statistic at its boundary value 0.
* The concentrated fixed-effects fits solve each unit's intercept exactly.
  The unit objective is a saw-tooth in the intercept (the support floor moves
  as the rate crosses integer multiples of `1 - alpha`, and the objective drops
  at every such breakpoint), so its supremum sits at an interior critical point
  of a tooth or at the left limit of a breakpoint; the inner search enumerates
  the breakpoints within a window of the envelope's peak (as wide as the
  widest tooth in the unit) and refines inside whole teeth by golden section;
  every evaluation searches each unit afresh, so the concentrated likelihood is
  a function of the parameters alone. It is continuous in the dispersion
  parameter but can jump in the slopes where a unit sits on its feasibility
  floor, so the outer optimizer is a deterministic multistart: BFGS with the
  analytic envelope gradient (the breakpoint a unit sits at is tracked in the
  gradient), a Nelder-Mead polish, and two perturbed restarts; `?cpb_fe`
  states the size of the alternative optima. The Katz family below `delta = 1` has a finite
  support too, but its likelihood only kinks at an integer ceiling (the
  entering point's mass grows from zero), so `gec_fe()` maximizes each unit's
  unimodal objective by golden section on an expanding bracket, with the
  feasibility floor exact (support `0..ceiling(mu/(1-delta))`). The compiled
  code asserts the
  unit index against the design, and the fits use only the estimation rows, so
  rows dropped for missing values no longer misalign the unit blocks.
* `cpb()` documents the discontinuity of the CPB log-likelihood and the
  deterministic multistart that maximizes it (BFGS from the Poisson solution
  at each `alpha.start`, Nelder-Mead polish, feasibility-repaired restarts);
  the mean-parameterized COM-Poisson tabulates its log-factorials.
* The CPB pmf is normalized by a recursion outward from the mode of its terms
  (consecutive terms differ by the factor `(n - k)/(k + 1) (1 - alpha)/alpha`,
  so the sum runs over the body of the distribution and stops where the terms
  fall below relative precision on either side), shared by the pooled
  likelihoods and the concentrated fixed-effects likelihood and gradient. A
  likelihood evaluation costs the width of the pmf rather than its ceiling: on
  counts in the hundreds one evaluation of the concentrated objective takes
  0.04 s instead of 2.6 s, and the pooled objective is 45 times cheaper. The
  R implementation of the pmf keeps the direct lgamma form as the independent
  check (agreement to 5e-10 across alpha from 0.01 to 0.999 and rates from
  0.05 to 2,000).
* The reported maximum of a pooled `cpb()` fit is the maximum of the fit's own
  profile in `alpha`: the profile is traced by continuation on both sides of the
  estimate (the slopes re-maximized at each step from the neighbouring solution
  and from the Poisson slopes, with a finer pass next to the estimate) until it
  has dropped four log-likelihood units, and a trace point above the
  multistart's value restarts the fit from there; the multistart itself runs
  from nine values of `alpha` and ends with an exact breakpoint search in
  `alpha` at the fitted slopes, since the teeth of the surface in `alpha` are
  thinner than any grid. `confint()` and `summary()`
  read the profile interval off the stored trace. On small samples the
  multistart alone could stop on a lower tooth with `alpha` off by 0.1; the
  remaining differences against a dense grid of starts are of the size of the
  jumps (below 0.005 log-likelihood units on 12 samples of 72 zero-truncated
  counts). Unit intercepts belong in `cpb_fe()`, which solves them exactly; a
  pooled fit with many dummy columns can stop short of it (0.14
  log-likelihood units on pscl's `prussian` with 14 corps dummies).
* The inner search of `cpb_fe()` locates the peak of a unit's smooth part to
  1e-6 (in 0.1.1's first candidate the twelve golden-section steps left an
  error of 1.6e-3, wider than a tooth at ceilings above a few hundred, so the
  wrong tooth could be searched) and searches the neighbouring tooth when the
  interior optimum sits at a tooth's end.
* The GEC support guard refuses instead of renormalizing: a rate whose recursion
  has not reached the tail by `max.support` is infeasible in the likelihood, a
  fit whose fitted support reaches the guard warns and records
  `$support_binding`, and `predict()`, `dgec()`, and `pgec()` return `NA` with
  a warning for such a rate (in 0.1.0 the pmf was silently renormalized over
  `0..max.support`, so a prediction far outside the fitted range was too small).
* The COM-Poisson support rule no longer clips the rate at 1e10 (with `mu =` at
  strong underdispersion the rate exceeds it and the pmf was wrong), and the
  d/p/q/r functions return `NA` for a non-finite rate.
* `.ud_nm_polish()` chains its restarts until a restart gains less than the
  tolerance, as documented (the first candidate stopped after one).
* Design matrices are rank-checked (an empty factor level or an exact linear
  combination is refused with the column named; a covariate constant within
  units is refused by the fixed-effects fits); a wrong-length vector `offset`,
  `weights`, or `cluster` is refused; `cluster` supplied without a bootstrap is
  announced; a cluster bootstrap with one cluster is refused and with fewer
  than thirty warned.
* `score()` raises `kmax` to the largest held-out count, and its predicted
  distributions carry the fit's offset (the in-sample log score equals
  `-logLik/n` for every class); the hurdle intensities are zero-truncated over
  the full support; `first_difference()` holds an offset at its mean and refuses
  a variable that enters an interaction; the Katz recursion runs at least to
  the observed count before its tail is dropped, so far-tail counts under an
  unbounded member get their small probability rather than an infeasibility
  sentinel; integer ceilings are feasible (floor with a tolerance) in the C++
  and the R pmf alike; the cluster-robust sandwich carries the G/(G-1) factor.
* The two-part models align both equations on the complete cases of every
  model variable; `zi_cpb()` records the user's call; `predict()` on `gec` and
  `gec_fe` uses the stored factor levels; the exported d/p/q/r functions
  return `numeric(0)` for zero-length input and refuse a vector dispersion
  parameter.
* Every estimator drops unused factor levels from its model frame, as `lm()`
  and `glm()` do, so data subset from a larger set no longer trip the
  rank-deficiency refusal on an empty level.
* A response whose counts are all zero still fits, with a warning that the rate
  runs to its lower bound and the dispersion parameter is not identified.
* `max.support`, `inner_it`, `alpha.start`, and `B` (with `se = "bootstrap"`)
  are checked on entry and refused by name before any fitting: a support below
  1 or beyond the integer range, a fractional or negative count, fewer than two
  bootstrap replicates, a start outside (0, 1).
* The support guard's warning comes once, from the fit the user asked for: the
  refits of a bootstrap, a bootstrap first difference, or the split-panel
  jackknife, the start values of the zero-inflated fits, the nest that
  `zi_test()` fits from a formula, and the hurdle arm of `compare_dispersion()`
  no longer repeat it, and a `gec()` or `gec_fe()` fit no longer issues it
  twice.

## Package

* `weights`, `cores`, and unknown-argument checks on the matched families
  (misspelled arguments are errors); `inst/CITATION`; a pkgdown site at
  https://bagozzib.github.io/underdisp/ and GitHub Actions running the full
  test suite (including the batteries skipped on CRAN) on five platforms.

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
* Fit-rank-exact degrees of freedom in the at-risk screen and the mixture
  comparisons, so zero-heavy designs are neither over- nor under-penalized.

## Data and integration

* Bundled `peacekeeping` dataset (UN peacekeeping contributions, a
  state-year roster with structural zeros).
* `broom` `tidy()`/`glance()` methods for all model classes (`augment()` for `cpb`), and
  `modelsummary` / `texreg` table support.
