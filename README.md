# underdisp

<!-- badges: start -->
<!-- badges: end -->

Tools for **detecting and modeling underdispersion** in count data — the case
where the conditional variance falls below the conditional mean, so counts
cluster more tightly around their expectation than a Poisson allows.

Underdispersion is common in bounded political-science counts (portfolios of
statuses that are filled and vacated over time) but poorly served by standard
software: the negative binomial cannot represent a variance below the mean and
collapses onto the Poisson. `underdisp` provides the missing pieces.

## What it does

- **Screening** — `ud_screen()` reports within-unit dispersion, a
  zero-truncated-Poisson *at-risk* benchmark (which separates genuine
  underdispersion from the artifact of conditioning on positive counts), and an
  over-conditioning guard, with calibrated or parametric-bootstrap thresholds.
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
- **Matched baselines** — `count_reg()`, `hurdle_count()`, and `zi_count()` fit
  Poisson, negative-binomial, and COM-Poisson analogues through the same
  interface, so `compare_models()` and `compare_dispersion()` adjudicate the
  whole family on one footing (information criteria plus proper scores, in
  sample, held out, or by cross-validation via `score()`/`cv_score()`).
- **Bias-corrected fixed effects** — `bias_correct = "jackknife"` removes the
  1/T incidental-parameters bias of the fixed-effects dispersion estimate
  (split-panel jackknife, Dhaene & Jochmans 2015), behind a validity gate that
  refuses the correction — with an informative warning — on panels that violate
  the method's time-homogeneity requirement.
- **Calibration and diagnostics** — hanging rootograms, PIT histograms, and
  `simulate()` methods for every model class, so any fit plugs into `DHARMa`'s
  simulated-residual diagnostics via `DHARMa::createDHARMa()`.
- **Quantities of interest** — predicted count distributions, the implied
  ceiling, incidence-rate ratios, and King-style first differences (including
  the exact extensive/intensive decomposition for two-part models), with
  bootstrap or profile-likelihood uncertainty; `broom`, `modelsummary`, and
  `texreg` support throughout.

## Installation

```r
# development version
# install.packages("remotes")
remotes::install_github("bagozzib/underdisp")

# from CRAN, once accepted
install.packages("underdisp")
```

## Quick start

```r
library(underdisp)

# an underdispersed count (var/mean ~ 0.5)
n <- 400; x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
d <- data.frame(y = rbinom(n, N, 0.5), x = x)

ud_screen(y ~ x, data = d)              # screen
fit <- cpb(y ~ x, data = d[d$y > 0, ])  # fit the zero-truncated CPB
summary(fit)
implied_ceiling(fit, newdata = data.frame(x = 0))
compare_dispersion(y ~ x, data = d)$table
```

The bundled `peacekeeping` panel demonstrates the package's central move — a
count that looks overdispersed in the pooled margin but is underdispersed
within countries at risk. See `vignette("underdisp")` for the full walk-through,
including the two-part models, the bias-corrected fixed effects, and the
DHARMa workflow.

## References

Dhaene, Geert, and Koen Jochmans. 2015. "Split-Panel Jackknife Estimation of
Fixed-Effect Models." *The Review of Economic Studies* 82(3): 991–1030.

King, Gary. 1989. "Variance Specification in Event Count Models." *American
Journal of Political Science* 33(3): 762–784.

Winkelmann, Rainer, Curtis S. Signorino, and Gary King. 1995. "A Correction for
an Underdispersed Event Count Probability Distribution." *Political Analysis* 5:
215–228.
