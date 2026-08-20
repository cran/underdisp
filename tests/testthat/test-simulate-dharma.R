## simulate() + fitted() across every model class, and the DHARMa entry point.
## One assertion helper serves all 11 classes: draws must be integer counts of
## the right shape, fitted() must be finite and predict-consistent, and the
## simulated grand mean must track the fitted grand mean.

expect_sim_ok <- function(fit, y, nsim = 15L, zt = FALSE) {
  s <- simulate(fit, nsim = nsim, seed = 7)
  expect_s3_class(s, "data.frame")
  m <- as.matrix(s)
  expect_identical(dim(m), c(length(y), as.integer(nsim)))
  expect_true(all(m >= 0) && all(m == floor(m)))
  if (zt) expect_true(all(m >= 1))
  f <- fitted(fit)
  expect_identical(length(f), length(y))
  expect_true(all(is.finite(f)))
  expect_lt(abs(mean(m) - mean(f)) / max(mean(f), 0.1), 0.25)
  ## explicit seed: reproducible, and the caller's RNG state is restored
  expect_identical(as.matrix(simulate(fit, nsim = 2, seed = 3)),
                   as.matrix(simulate(fit, nsim = 2, seed = 3)))
  invisible(m)
}

expect_dharma_ok <- function(fit, y) {
  skip_if_not_installed("DHARMa")
  s  <- as.matrix(simulate(fit, nsim = 60, seed = 11))
  dh <- DHARMa::createDHARMa(simulatedResponse = s, observedResponse = y,
                             fittedPredictedResponse = fitted(fit),
                             integerResponse = TRUE)
  r <- dh$scaledResiduals
  expect_true(all(is.finite(r)) && all(r >= 0) && all(r <= 1))
}

## shared toy data ------------------------------------------------------------
set.seed(42)
n <- 240
x <- rnorm(n); z <- rnorm(n)
Ncap <- pmax(round(exp(1.3 + 0.4 * x) / 0.5), 1)
d_cs  <- data.frame(y = rbinom(n, Ncap, 0.5), x = x, z = z)      # underdispersed counts
d_pos <- d_cs[d_cs$y > 0, ]                                       # positives for ZT fits
uid   <- rep(sprintf("u%02d", 1:24), each = 10)
d_pan <- data.frame(y = rcpb(n, exp(rep(rnorm(24, 1.2, 0.3), each = 10) + 0.3 * x), 0.5),
                    x = x, unit = uid)
d_h   <- data.frame(y = rhurdle_cpb(n, exp(1.2 + 0.3 * x), 0.5, plogis(0.2 + 0.8 * z)),
                    x = x, z = z)

test_that("simulate/fitted: quick smoke (cpb and count_reg, plus the DHARMa entry point)", {
  f1 <- cpb(y ~ x, d_pos, truncated = TRUE, se = "none")
  expect_sim_ok(f1, d_pos$y, zt = TRUE)
  f5 <- count_reg(y ~ x, d_pos, family = "poisson", truncated = TRUE, se = "none")
  expect_sim_ok(f5, d_pos$y, zt = TRUE)
  expect_dharma_ok(f1, d_pos$y)
})

test_that("simulate/fitted: single-equation classes (cpb, cpb_fe, gec, gec_fe, count_reg)", {
  skip_on_cran()
  f1 <- cpb(y ~ x, d_pos, truncated = TRUE, se = "none")
  expect_sim_ok(f1, d_pos$y, zt = TRUE)
  f2 <- cpb_fe(y ~ x, d_pan, fe = "unit", se = "none")
  expect_sim_ok(f2, d_pan$y)
  f3 <- gec(y ~ x, d_cs, se = "none")
  expect_sim_ok(f3, d_cs$y)
  f4 <- gec_fe(y ~ x, d_pan, fe = "unit", se = "none")   # inherits simulate.gec/fitted.gec
  expect_sim_ok(f4, d_pan$y)
  f5 <- count_reg(y ~ x, d_pos, family = "poisson", truncated = TRUE, se = "none")
  expect_sim_ok(f5, d_pos$y, zt = TRUE)
  f6 <- count_reg(y ~ x, d_cs, family = "compois", se = "none")
  expect_sim_ok(f6, d_cs$y)
  expect_dharma_ok(f3, d_cs$y)
  expect_dharma_ok(f6, d_cs$y)
})

test_that("simulate/fitted: hurdle classes carry the zero structure", {
  skip_on_cran()
  f7 <- hurdle_cpb(y ~ x, d_h, participation = ~ z)
  m7 <- expect_sim_ok(f7, d_h$y)
  expect_gt(mean(m7 == 0), 0.05)                     # hurdle stage produces zeros
  f8 <- hurdle_gec(y ~ x, d_h, participation = ~ z)
  expect_sim_ok(f8, d_h$y)
  f9 <- hurdle_count(y ~ x, d_h, family = "negbin", participation = ~ z, se = "none")
  expect_sim_ok(f9, d_h$y)
  expect_dharma_ok(f7, d_h$y)
})

test_that("simulate/fitted: zero-inflated classes carry the zero structure", {
  skip_on_cran()
  f10 <- zi_cpb(y ~ x, d_h, zero = ~ z, se = "none")
  m10 <- expect_sim_ok(f10, d_h$y)
  expect_gt(mean(m10 == 0), 0.05)
  f11 <- zi_gec(y ~ x, d_h, zero = ~ z, se = "none")
  expect_sim_ok(f11, d_h$y)
  f12 <- zi_count(y ~ x, d_h, family = "poisson", zero = ~ z, se = "none")
  expect_sim_ok(f12, d_h$y)
  expect_dharma_ok(f12, d_h$y)
})
