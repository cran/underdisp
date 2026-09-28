# underdisp

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/underdisp)](https://CRAN.R-project.org/package=underdisp)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/underdisp)](https://CRAN.R-project.org/package=underdisp)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/bagozzib/underdisp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bagozzib/underdisp/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bagozzib/underdisp/actions/workflows/pkgdown.yaml/badge.svg)](https://bagozzib.github.io/underdisp/)
<!-- badges: end -->

Tools for **detecting and modeling underdispersion** in count data: the case
where the conditional variance falls below the conditional mean, so counts
cluster more tightly around their expectation than a Poisson allows.

Underdispersion is common in bounded counts (portfolios of statuses that are
filled and vacated over time, events with regular timing, quotas) but poorly
served by standard software: the negative binomial cannot represent a variance
below the mean and collapses onto the Poisson. `underdisp` provides the missing
pieces, from a screen that tells you whether an outcome is underdispersed to a
family of estimators, tests, and diagnostics that treat every model on the same
footing.

## What it does

- **Screening** — `ud_screen()` reports within-unit dispersion, a
  zero-truncated-Poisson *at-risk* benchmark (which separates genuine
  underdispersion from the artifact of conditioning on positive counts), and an
  over-conditioning guard, with calibrated or parametric-bootstrap thresholds.
- **Tests and profiles** — `dispersion_test()` tests any fitted family's
  dispersion parameter against the Poisson on the fitted design (fixed effects
  included), with boundary-aware p-values; `dispersion_profile()` plots the
  conditional variance-to-mean ratio against the fitted mean next to the curve
  each family implies, so the mechanism (a hard ceiling, regular event timing,
  a soft tail) can be read off the data.
- **The hard-ceiling family** — `cpb()` fits King's continuous parameter
  binomial and its zero-truncated variant, with an interpretable
  observation-specific bound; `cpb_fe()` absorbs high-dimensional unit fixed
  effects by a concentrated likelihood that scales to thousands of units.
- **The free-dispersion family** — `gec()` fits the generalized event count
  (Katz) model, whose single dispersion parameter is estimated freely and spans
  under-, equi-, and overdispersion; `gec_fe()` adds concentrated fixed effects.
- **Excess zeros** — `hurdle_cpb()`/`zi_cpb()` and `hurdle_gec()`/`zi_gec()`
  pair the underdispersed intensities with participation or structural-zero
  processes; `zi_test()` runs the boundary-corrected zero-inflation test.
- **Matched families** — `count_reg()`, `hurdle_count()`, and `zi_count()` fit
  the Poisson, negative binomial, COM-Poisson (rate- or mean-parameterized),
  generalized Poisson, gamma-count, and double Poisson through the same
  interface, with fixed effects, offsets, frequency weights, and analytic,
  robust, or cluster-robust standard errors, so `compare_models()` and
  `compare_dispersion()` adjudicate the whole family on one footing
  (information criteria plus proper scores, in sample, held out, or by
  cross-validation via `score()`/`cv_score()`).
- **Bias-corrected fixed effects** — `bias_correct = "jackknife"` removes the
  1/T incidental-parameters bias of the fixed-effects dispersion estimate
  (split-panel jackknife, Dhaene & Jochmans 2015), behind a validity gate that
  refuses the correction, with an informative warning, on panels that violate
  the method's time-homogeneity requirement.
- **Calibration and diagnostics** — hanging rootograms, PIT histograms, and
  `simulate()` methods for every model class, so any fit plugs into `DHARMa`'s
  simulated-residual diagnostics via `DHARMa::createDHARMa()`.
- **Quantities of interest** — predicted count distributions, the implied
  ceiling, incidence-rate ratios, and first differences (including the exact
  extensive/intensive decomposition for two-part models), with bootstrap,
  delta-method, or profile-likelihood uncertainty; `confint()` for every
  class; parallel bootstraps via `cores =`; `broom`, `modelsummary`, and
  `texreg` support throughout.

## Installation

```r
# from CRAN
install.packages("underdisp")

# development version
# install.packages("remotes")
remotes::install_github("bagozzib/underdisp")
```

## Quick start

```r
library(underdisp)

# an underdispersed count (var/mean ~ 0.5)
n <- 400; x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
d <- data.frame(y = rbinom(n, N, 0.5), x = x)

ud_screen(y ~ x, data = d)                     # screen
fit <- cpb(y ~ x, data = d, truncated = FALSE, se = "none")
summary(fit)
dispersion_test(fit, cores = 2)                # boundary null: parametric-bootstrap p-value
implied_ceiling(fit, newdata = data.frame(x = 0))
compare_dispersion(y ~ x, data = d)$table
dispersion_profile(cpb = fit,                  # which mechanism fits?
                   gammacount = count_reg(y ~ x, d, family = "gammacount"),
                   negbin = count_reg(y ~ x, d, family = "negbin"))
```

The bundled `peacekeeping` panel demonstrates the package's central move: a
count that looks overdispersed in the pooled margin but is underdispersed
within countries at risk. See `vignette("underdisp")` for the full walk-through,
including the two-part models, the bias-corrected fixed effects, the wider
family of underdispersed distributions, and the DHARMa workflow. The package
website is at <https://bagozzib.github.io/underdisp/>.

## References

Dhaene, Geert, and Koen Jochmans. 2015. "Split-Panel Jackknife Estimation of
Fixed-Effect Models." *The Review of Economic Studies* 82(3): 991–1030.

Efron, Bradley. 1986. "Double Exponential Families and Their Use in Generalized
Linear Regression." *Journal of the American Statistical Association* 81(395):
709–721.

Huang, Alan. 2017. "Mean-Parametrized Conway–Maxwell–Poisson Regression Models
for Dispersed Counts." *Statistical Modelling* 17(6): 359–380.

King, Gary. 1989. "Variance Specification in Event Count Models." *American
Journal of Political Science* 33(3): 762–784.

Winkelmann, Rainer. 1995. "Duration Dependence and Dispersion in Count-Data
Models." *Journal of Business & Economic Statistics* 13(4): 467–474.

Winkelmann, Rainer, Curtis S. Signorino, and Gary King. 1995. "A Correction for
an Underdispersed Event Count Probability Distribution." *Political Analysis* 5:
215–228.
