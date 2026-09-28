## Frequency weights: a weight of w equals w copies of the row, in every estimator.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

set.seed(7); n <- 300; x <- rnorm(n); z <- rnorm(n)
d0 <- data.frame(y = rcpb(n, exp(1.2 + 0.5 * x), 0.5), x = x, z = z, u = factor(sample(1:15, n, TRUE)))
w  <- sample(1:3, n, TRUE); ddup <- d0[rep(seq_len(n), w), ]
same_fit <- function(f1, f2, tol_ll = 1e-3, tol_b = 5e-3) {
  expect_equal(f1$loglik, f2$loglik, tolerance = tol_ll)
  expect_equal(unlist(coef(f1)), unlist(coef(f2)), tolerance = tol_b)
}

test_that("weights equal duplicated rows for the single-equation estimators", {
  same_fit(cpb(y ~ x, d0, truncated = FALSE, se = "none", weights = w), cpb(y ~ x, ddup, truncated = FALSE, se = "none"), 0.05)
  same_fit(gec(y ~ x, d0, se = "none", weights = w), gec(y ~ x, ddup, se = "none"), 0.05)
  for (fam in c("poisson", "negbin", "compois"))
    same_fit(count_reg(y ~ x, d0, family = fam, weights = w), count_reg(y ~ x, ddup, family = fam))
  cw <- count_reg(y ~ x, d0, family = "poisson", weights = w); cd <- count_reg(y ~ x, ddup, family = "poisson")
  expect_equal(sqrt(diag(vcov(cw))), sqrt(diag(vcov(cd))), tolerance = 1e-3)   # information as duplicated rows
  expect_equal(nobs(cw), sum(w)); expect_equal(nobs(cd), sum(w))
  expect_equal(as.numeric(logLik(cw)), as.numeric(logLik(cd)), tolerance = 1e-6)
  cr <- count_reg(y ~ x, data = transform(d0, wcol = w), family = "poisson", weights = "wcol")   # column-name form
  expect_equal(cr$loglik, cw$loglik)
})

test_that("weights equal duplicated rows for the fixed-effects and two-part estimators", {
  skip_on_cran()
  same_fit(cpb_fe(y ~ x, d0, fe = "u", weights = w), cpb_fe(y ~ x, ddup, fe = "u"), 0.1, 1e-2)
  same_fit(gec_fe(y ~ x, d0, fe = "u", weights = w), gec_fe(y ~ x, ddup, fe = "u"), 0.1, 1e-2)
  ## the two-part models need real structural zeros, or their zero part is unidentified
  set.seed(8); pz <- plogis(-0.8 + 0.6 * z); dz <- d0; dz$y <- ifelse(runif(n) < pz, 0L, d0$y)
  dzdup <- dz[rep(seq_len(n), w), ]
  same_fit(hurdle_cpb(y ~ x, dz, participation = ~ z, weights = w), hurdle_cpb(y ~ x, dzdup, participation = ~ z), 0.1, 1e-2)
  same_fit(hurdle_count(y ~ x, dz, family = "negbin", participation = ~ z, weights = w),
           hurdle_count(y ~ x, dzdup, family = "negbin", participation = ~ z))
  same_fit(zi_count(y ~ x, dz, family = "poisson", zero = ~ z, weights = w), zi_count(y ~ x, dzdup, family = "poisson", zero = ~ z))
  ## the CPB mixtures are maximized on a discontinuous surface (see ?cpb): the two
  ## runs may settle on different teeth, so the comparison is at the tooth scale
  same_fit(zi_cpb(y ~ x, dz, zero = ~ z, weights = w), zi_cpb(y ~ x, dzdup, zero = ~ z), 0.5, 0.1)
  same_fit(zi_gec(y ~ x, dz, zero = ~ z, weights = w), zi_gec(y ~ x, dzdup, zero = ~ z), 0.5, 0.1)
})

test_that("bad weights are refused with the argument named", {
  expect_error(cpb(y ~ x, d0, truncated = FALSE, se = "none", weights = w[-1]), "weights")
  expect_error(count_reg(y ~ x, d0, weights = -w), "non-negative")
  expect_error(gec(y ~ x, d0, se = "none", weights = "nope"), "weights")
  expect_error(zi_cpb(y ~ x, d0, zero = ~ z, weights = w, method = "em"), "weights")
})
