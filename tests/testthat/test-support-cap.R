## The CPB and GEC ceilings are lambda / (1 - alpha) and mu / (1 - delta). Both diverge as the dispersion
## parameter approaches its equidispersed limit, and the likelihood guards only bar alpha within 1e-12 of 1 --
## which still leaves a ceiling of 1e12, far outside int's range. Converting it was undefined behaviour, and
## the hardware disagreed about the result (x86-64 gives INT_MIN, arm64 saturates to INT_MAX). CRAN's M1
## sanitizer build reported it on 2026-09-28 at cpb_pmf.h:30, cpb_fe.cpp:106 and cpb_fe.cpp:164.
##
## These tests pin the behaviour the cap guarantees: a ceiling past max_support is infeasible, on every
## platform, and fits whose dispersion is driven to the limit still return rather than acting on a wrapped
## integer. They exercise the C++ kernel directly, because the d/p/q functions compute the pmf in R and so
## would not reach the patched code at all.

test_that("the likelihood kernel rejects every ceiling past max_support, however large", {
  X <- matrix(1, 40, 1); Y <- rep(1L, 40); off <- rep(0, 40)
  nll <- function(alpha) underdisp:::cpb_nll_cpp(c(0, log(alpha / (1 - alpha))), X, Y, off, 500L, FALSE)

  ## feasible ceilings give a finite objective
  expect_equal(nll(0.5), 27.725887222397784, tolerance = 1e-12)
  expect_equal(nll(0.9), 37.929785636817449, tolerance = 1e-12)

  ## ceilings beyond max_support are infeasible and return the sentinel, including the ones that overflow
  ## int: 1e9 is merely past max_support, while 5e11 and 1e13 are past int's range entirely.
  for (a in c(1 - 1e-9, 1 - 2e-12, 1 - 1e-12, 1 - 1e-13))
    expect_identical(nll(a), 1e10)
})

test_that("a ceiling that is an exact integer keeps its own count feasible", {
  ## lambda / (1 - alpha) = 10 exactly: 10 is inside the support, 11 is not.
  expect_true(is.finite(dcpb(10L, 5, 0.5, log = TRUE)))
  expect_identical(dcpb(11L, 5, 0.5, log = TRUE), -Inf)
})

test_that("fits whose dispersion is driven to the equidispersed limit still return", {
  ## Overdispersed data give alpha nowhere to go but its upper limit, which is where the ceiling diverges.
  ## The warning about the support guard is the expected outcome here, not the thing under test.
  set.seed(4); n <- 200; x <- rnorm(n)
  y <- rnbinom(n, mu = exp(1.2 + 0.4 * x), size = 1.5)
  f <- suppressWarnings(cpb(y ~ x, data.frame(y = y, x = x), truncated = FALSE, se = "none"))
  expect_true(is.finite(f$loglik))
  expect_true(f$alpha > 0 && f$alpha <= 1)
  expect_true(all(is.finite(coef(f))))

  ## the same on the fixed-effects path, which carries four of the six casts, with a Poisson panel
  set.seed(5); NU <- 25; TT <- 8; id <- rep(seq_len(NU), each = TT)
  xp <- rnorm(NU * TT); yp <- rpois(NU * TT, exp(1 + 0.3 * xp + rep(rnorm(NU, 0, 0.4), each = TT)))
  pan <- data.frame(y = yp, x = xp, id = id)
  ffe <- suppressWarnings(cpb_fe(y ~ x, pan, fe = "id", se = "none"))
  expect_true(is.finite(ffe$loglik))
  expect_true(all(is.finite(ffe$fe)))

  ## and on the Katz side, whose tooth scan carries the sixth cast
  gfe <- suppressWarnings(gec_fe(y ~ x, pan, fe = "id", se = "none"))
  expect_true(is.finite(gfe$loglik))
})
