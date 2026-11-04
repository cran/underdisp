## The separation check of the zero-inflation equation (R/zi_separation.R): the profile of each covariate
## coefficient, the rest of the zero equation re-maximized with the count component held at its estimate.

zi_sep_designs <- function() {
  set.seed(11); n <- 400; x <- rnorm(n); z <- rnorm(n); zb <- rbinom(n, 1, 0.3)
  list(
    clean     = data.frame(y = ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x))), x, z, zb),
    ## every observation with zb = 1 is a zero: the classic quasi-separation
    classic   = data.frame(y = ifelse(zb == 1, 0L, ifelse(rbinom(n, 1, 0.2) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))), x, z, zb),
    ## zeros sorted by a threshold in z, the zeros below it sampling zeros: separation of the mixture's own kind,
    ## with no covariate splitting the zeros from the positives
    threshold = data.frame(y = ifelse(z > 0.8, 0L, rpois(n, exp(0.2 + 0.3 * x))), x, z, zb),
    ## zb = 1 makes a structural zero likely (0.85) but leaves positives there: finite
    finite    = { set.seed(12); data.frame(y = ifelse(rbinom(n, 1, ifelse(zb == 1, 0.85, 0.2)) == 1, 0L,
                                                     rpois(n, exp(1 + 0.3 * x))), x, z, zb) })
}

test_that("the check flags the coefficients with no finite estimate and only those", {
  skip_on_cran()
  d <- zi_sep_designs()
  sep <- function(dd, zero) {
    f <- suppressWarnings(zi_count(y ~ x, dd, family = "poisson", zero = zero, se = "none"))
    underdisp:::.zi_separation(f$Zz, f$zero.coefficients, f$y, dpois(0, f$lambda_full))$diverging
  }
  expect_identical(sep(d$clean, ~ z), character(0))
  expect_identical(sep(d$classic, ~ zb + z), "zb")
  expect_identical(sep(d$threshold, ~ z), "z")
  expect_identical(sep(d$finite, ~ zb + z), character(0))
})

test_that("short_of_max measures from the fit, so a fit off the maximum shows its shortfall", {
  d <- zi_sep_designs()
  f <- suppressWarnings(zi_count(y ~ x, d$clean, family = "poisson", zero = ~ z, se = "none"))
  f0 <- dpois(0, f$lambda_full)
  expect_lt(underdisp:::.zi_separation(f$Zz, f$zero.coefficients, f$y, f0)$short_of_max, 0.01)
  ## the intercept moved off the maximum: the profile of z at its estimate re-maximizes the intercept and returns
  ## to the fit's maximum, so the shortfall found is the whole gap
  g <- f$zero.coefficients; g["(Intercept)"] <- g["(Intercept)"] + 0.4
  lz <- function(g) { p <- plogis(as.numeric(f$Zz %*% g)); sum(ifelse(f$y == 0, log(p + (1 - p) * f0), log1p(-p))) }
  gap <- lz(f$zero.coefficients) - lz(g)
  expect_gt(gap, 0.5)
  expect_equal(underdisp:::.zi_separation(f$Zz, g, f$y, f0)$short_of_max, gap, tolerance = 1e-3)
})

test_that("the restart keeps a better point, checks again, and stops when it no longer pays", {
  ## a stand-in fit whose log-likelihood the check can raise twice, by 1 and then by 0.5
  gains <- c(1, 0.5, 0)
  check <- function(f) list(short_of_max = gains[f$k + 1L], best = if (gains[f$k + 1L] > 0) f$k + 1L)
  refit <- function(f, g) list(k = g, ll = f$ll + gains[f$k + 1L])
  r <- underdisp:::.zi_restart(list(k = 0L, ll = -10), check, refit, function(f) f$ll)
  expect_identical(r$separation$restarts, 2L)
  expect_equal(r$separation$restart_gain, 1.5)
  expect_equal(r$fit$ll, -8.5)
  ## a restart that loses is not kept
  r2 <- underdisp:::.zi_restart(list(k = 0L, ll = -10), check, function(f, g) list(k = g, ll = f$ll - 1), function(f) f$ll)
  expect_identical(r2$separation$restarts, 0L)
  expect_equal(r2$fit$ll, -10)
})

test_that("zi_cpb's restart path re-polishes from a moved zero equation back to the maximum", {
  set.seed(1); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
  d <- data.frame(y, x, z)
  h <- zi_cpb(y ~ x, d, zero = ~ z)
  expect_identical(h$separation$restarts, 0L)
  g <- h$zero_coef; g["z"] <- g["z"] + 0.5
  r <- underdisp:::.zi_cpb_fit(y ~ x, d, ~ z, support = "ceiling",
                               restart = list(b = h$coefficients, alpha = h$alpha, g = g, iterations = h$iterations))
  expect_equal(r$loglik, h$loglik, tolerance = 1e-6)
  expect_equal(r$zero_coef, h$zero_coef, tolerance = 1e-3)
})

test_that("zi_count() warns, stores the check, and prints a note", {
  d <- zi_sep_designs()
  expect_warning(f <- zi_count(y ~ x, d$classic, family = "poisson", zero = ~ zb + z, se = "none"), "separated or nearly so")
  expect_identical(f$separation$diverging, "zb")
  expect_output(print(f), "zero-inflation equation is separated")
  expect_no_warning(g <- zi_count(y ~ x, d$finite, family = "poisson", zero = ~ zb + z, se = "none"))
  expect_length(g$separation$diverging, 0L)
})

test_that("a fit to data with no inflation is near the boundary, not separated", {
  ## no structural zeros: the zero equation's coefficients are unidentified because the fit sits at the
  ## no-inflation boundary, which every profile can reach, so nothing is flagged and the fit says why
  set.seed(21); n <- 300; x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(y = rpois(n, exp(0.8 + 0.4 * x)), x, z)
  expect_no_warning(f <- zi_count(y ~ x, d, family = "poisson", zero = ~ z, se = "none"))
  expect_length(f$separation$diverging, 0L)
  expect_true(f$separation$no_inflation)
  expect_output(print(f), "near the no-inflation boundary")
  ## and an inflated design is not reported as near it
  g <- suppressWarnings(zi_count(y ~ x, zi_sep_designs()$clean, family = "poisson", zero = ~ z, se = "none"))
  expect_false(g$separation$no_inflation)
})

test_that("a unit with no positive count is a boundary of its own effect, not separation", {
  set.seed(5); nu <- 30; tt <- 8; u <- rep(seq_len(nu), each = tt); x <- rnorm(nu * tt); z <- rnorm(nu * tt)
  y <- ifelse(rbinom(nu * tt, 1, plogis(-0.8 + 0.6 * z)) == 1, 0L, rpois(nu * tt, exp(0.8 + 0.3 * x)))
  y[u <= 5] <- 0L                                         # five units with no positive count
  d <- data.frame(y, x, z, u = factor(u))
  expect_no_warning(f <- zi_count(y ~ x, d, family = "poisson", zero = ~ z, zero_fe = "u", se = "none"))
  expect_length(f$separation$diverging, 0L)
})

test_that("zi_cpb() and zi_gec() run the same check", {
  skip_on_cran()                                          # each mixture fit takes a few seconds
  d <- zi_sep_designs()
  expect_warning(f <- zi_cpb(y ~ x, d$classic, zero = ~ zb + z), "zi_cpb: the zero-inflation equation is separated")
  expect_true("zb" %in% f$separation$diverging)
  expect_warning(g <- zi_gec(y ~ x, d$classic, zero = ~ zb + z), "zi_gec: the zero-inflation equation is separated")
  expect_true("zb" %in% g$separation$diverging)
  set.seed(1); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
  h <- zi_cpb(y ~ x, data.frame(y, x, z), zero = ~ z)
  expect_length(h$separation$diverging, 0L)
})
