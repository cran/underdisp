## King's CPB (support = "ceiling", the default from 0.2.0) is the GEC at delta = alpha < 1 (Winkelmann,
## Signorino and King 1995); the floor support is the definition of versions before 0.2.0. These tests pin
## the identity, the edge of the support, the Poisson limit, and the stored setting.

test_that("King's CPB pmf is the GEC's at delta = alpha", {
  for (a in c(0.1, 0.5, 0.9, 0.999)) for (l in c(0.2, 1, 3, 7.3, 60)) {
    k <- 0:min(ceiling(l / (1 - a)) + 2, ceiling(l + 60 * sqrt(l + 1) + 60))
    expect_equal(dcpb(k, l, a), dgec(k, l, a, max.support = 100000L), tolerance = 1e-12)
    expect_equal(dcpb(k[-1], l, a, truncated = TRUE),
                 dgec(k[-1], l, a, max.support = 100000L) / (1 - dgec(0, l, a, max.support = 100000L)),
                 tolerance = 1e-12)
  }
})

test_that("the C++ log pmf on King's support matches the R pmf", {
  for (a in c(0.2, 0.6, 0.95)) for (l in c(0.4, 2.6, 12)) {
    k <- which(dcpb(0:80, l, a) > 1e-12) - 1L
    m <- underdisp:::cpb_lp0_cpp(c(log(l), qlogis(a)), matrix(1, length(k), 1), as.integer(k), rep(0, length(k)),
                                 500L, FALSE, FALSE)
    expect_equal(m[, 1], log(dcpb(k, l, a)), tolerance = 1e-10)
    expect_equal(m[, 2], rep(dcpb(0, l, a), length(k)), tolerance = 1e-12)
  }
})

test_that("the two supports differ by the fractional top point and coincide at integer ceilings", {
  expect_equal(dcpb(0:7, 3, 0.5), dcpb(0:7, 3, 0.5, support = "floor"), tolerance = 1e-14)   # n = 6
  expect_gt(dcpb(6, 2.6, 0.5), 0)                                                           # n = 5.2
  expect_equal(dcpb(6, 2.6, 0.5, support = "floor"), 0)
  expect_equal(dcpb(7, 2.6, 0.5), 0)
  ## the realized mean sits at or above lambda on King's support (Winkelmann, Signorino and King 1995,
  ## Proposition 1) and at or below it on the floor support
  mn <- function(l, a, s) sum((0:200) * dcpb(0:200, l, a, support = s))
  for (l in c(0.7, 2.6, 4.1)) {
    expect_gte(mn(l, 0.5, "ceiling"), l - 1e-12)
    expect_lte(mn(l, 0.5, "floor"), l + 1e-12)
  }
})

test_that("a count may exceed lambda/(1-alpha) by less than one on King's support, and no more", {
  X <- matrix(1, 1, 1); par <- c(log(2.6), qlogis(0.5))                   # n = 5.2
  expect_lt(underdisp:::cpb_nll_cpp(par, X, 6L, 0, 500L, FALSE, FALSE), 1e9)
  expect_equal(underdisp:::cpb_nll_cpp(par, X, 7L, 0, 500L, FALSE, FALSE), 1e10)
  expect_equal(underdisp:::cpb_nll_cpp(par, X, 6L, 0, 500L, FALSE, TRUE), 1e10)   # floor support: 6 > 5.2
})

test_that("King's zero-truncated CPB exists for every positive rate", {
  expect_equal(sum(dcpb(1:5, 0.3, 0.5, truncated = TRUE)), 1, tolerance = 1e-12)
  set.seed(1)
  expect_silent(rcpb(50, 0.3, 0.5, truncated = TRUE))
  expect_warning(rcpb(50, 0.3, 0.5, truncated = TRUE, support = "floor"), "below 1")
})

test_that("King's derivatives match central differences of the log pmf", {
  set.seed(2)
  for (r in 1:60) {
    l <- exp(runif(1, -1, 3)); a <- runif(1, 0.05, 0.95); tr <- r %% 2 == 0
    nn <- l / (1 - a); if (abs(nn - round(nn)) < 1e-3) next
    y <- sample(if (tr) 1:ceiling(nn) else 0:ceiling(nn), 1)
    lp <- function(eta, al) log(dcpb(y, exp(eta), al, truncated = tr))
    h <- 1e-6
    an <- underdisp:::cpb_dlogpmf_cpp(as.integer(y), l, a, 500L, tr, FALSE)
    expect_equal(an[1, 1], (lp(log(l) + h, a) - lp(log(l) - h, a)) / (2 * h), tolerance = 1e-5)
    expect_equal(an[1, 2], (lp(log(l), a + h) - lp(log(l), a - h)) / (2 * h), tolerance = 1e-5)
  }
})

test_that("cpb() on King's support is gec() restricted to delta < 1", {
  set.seed(31); n <- 400; x <- rnorm(n)
  d <- data.frame(y = as.integer(rgec(n, exp(1 + 0.5 * x), 0.4)), x = x)
  fc <- cpb(y ~ x, d, truncated = FALSE); fg <- gec(y ~ x, d)
  expect_identical(fc$support, "ceiling")
  expect_equal(fc$loglik, fg$loglik, tolerance = 1e-6)
  expect_equal(fc$alpha, fg$delta, tolerance = 1e-3)
  expect_equal(unname(coef(fc)), unname(coef(fg)), tolerance = 1e-3)
  dp <- d[d$y > 0, ]
  expect_equal(cpb(y ~ x, dp)$loglik, gec(y ~ x, dp, truncated = TRUE)$loglik, tolerance = 1e-6)
})

test_that("cpb_fe() on King's support matches gec_fe(), and its zero-truncated fit is a maximum", {
  set.seed(4); U <- 25; Tn <- 8
  pd <- do.call(rbind, lapply(1:U, function(u) { x <- rnorm(Tn)
    data.frame(unit = u, x = x, y = as.integer(rgec(Tn, exp(rnorm(1, 1, 0.4) + 0.5 * x), 0.4))) }))
  f1 <- cpb_fe(y ~ x, pd, fe = "unit"); f2 <- gec_fe(y ~ x, pd, fe = "unit")
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-5)
  expect_equal(f1$alpha, f2$delta, tolerance = 1e-3)
  pp <- pd[pd$y > 0, ]; pp <- pp[order(pp$unit), ]
  ft <- cpb_fe(y ~ x, pp, fe = "unit", truncated = TRUE)
  lam <- exp(ft$fe[as.character(pp$unit)] + coef(ft)[1] * pp$x)
  expect_equal(sum(log(dcpb(pp$y, lam, ft$alpha, truncated = TRUE))), ft$loglik, tolerance = 1e-8)
  nll <- function(par) { u <- par[seq_len(U)]; lam <- exp(u[pp$unit] + par[U + 1] * pp$x)
    v <- -sum(log(dcpb(pp$y, lam, plogis(par[U + 2]), truncated = TRUE))); if (is.finite(v)) v else 1e10 }
  st <- c(unname(ft$fe[as.character(1:U)]), coef(ft)[1], qlogis(ft$alpha))
  op <- optim(st, nll, method = "BFGS", control = list(maxit = 200, reltol = 1e-12))
  expect_lt(-op$value - ft$loglik, 1e-6)
})

test_that("King's CPB approaches the Poisson as alpha -> 1, without the guard's cliff", {
  expect_equal(dcpb(0:30, 4, 1 - 1e-9), dpois(0:30, 4), tolerance = 1e-6)
  set.seed(32); d <- data.frame(y = stats::rnbinom(300, mu = 3, size = 4))     # overdispersed
  expect_warning(f <- cpb(y ~ 1, d, truncated = FALSE), "equidispersed limit")
  expect_gt(f$alpha, 0.999)
  expect_false(f$support_binding)
  pois <- sum(stats::dpois(d$y, mean(d$y), log = TRUE))
  expect_lt(abs(f$loglik - pois), 0.05)
})

test_that("the C++ kernel stays a probability near the Poisson limit and outside the pmf's body", {
  lp0 <- function(y, l, a) underdisp:::cpb_lp0_cpp(c(log(l), qlogis(a)), matrix(1, length(y), 1), as.integer(y),
                                                   rep(0, length(y)), 100000L, FALSE, FALSE)
  ## near the limit the trials parameter lambda / (1 - alpha) reaches 1e13; P(0) of a large rate lies outside the
  ## summed body and is taken in closed form (it once came out above 1, and a zero-inflated fit then reported a
  ## positive log-likelihood)
  for (l in c(3, 60, 200)) for (a in 1 - 10^-(6:11)) {
    ky <- c(0, 1, round(l)); m <- lp0(ky, l, a)
    expect_true(all(m[, 1] <= 0) && all(m[, 2] <= 1))
    expect_equal(m[, 1], stats::dpois(ky, l, log = TRUE), tolerance = 1e-4)
  }
  ## within 1e-12 of one: infeasible, as in the likelihoods. The kernel's sentinel is the compiler's -1e300, which
  ## R's own parse of the literal can miss by an ulp where long double is double (macOS arm64), so the test reads
  ## it the way the package does, as anything below -1e299
  expect_lt(lp0(5, 200, 1 - 1e-15)[1, 1], -1e299)
  ## counts far below a narrow body (a large rate, a tiny alpha): closed form, against a direct lgamma evaluation
  l <- 2000; a <- 0.001; n <- l / (1 - a); k <- 0:ceiling(n - 1e-9)
  lt <- lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + k * log(1 - a) + (n - k) * log(a)
  ref <- lt - (max(lt) + log(sum(exp(lt - max(lt)))))
  expect_equal(lp0(c(1, 100, 1000), l, a)[, 1], ref[c(1, 100, 1000) + 1], tolerance = 1e-10)
  ## and without walking to them: a rate of 1e8 with alpha 1e-9 is a body of a few points 1e8 steps from the count
  expect_lt(system.time(for (i in 1:50) lp0(1, 1e8, 1e-9))[["elapsed"]], 1)
})

test_that("fits store the support, and objects saved before 0.2.0 read as floor", {
  set.seed(33); d <- data.frame(y = rcpb(200, 3.3, 0.5, truncated = TRUE))
  f <- cpb(y ~ 1, d)
  expect_identical(f$support, "ceiling")
  ff <- cpb(y ~ 1, d, support = "floor")
  expect_identical(ff$support, "floor")
  old <- ff; old$support <- NULL
  expect_identical(underdisp:::.cpb_support(old), "floor")
  expect_equal(predict(old, type = "response"), predict(ff, type = "response"))   # the methods read the support
  kold <- f; kold$support <- NULL                                   # a King's fit stripped of it reads as floor
  expect_false(isTRUE(all.equal(predict(kold, type = "response"), predict(f, type = "response"))))
  expect_error(cpb(y ~ 1, d, support = "round"))
})

test_that("a King's zi_cpb predicts with the exact count mean; the floor support keeps the rate", {
  skip_on_cran()
  set.seed(34); n <- 400; x <- rnorm(n)
  d <- data.frame(y = rzicpb(n, exp(0.9 + 0.4 * x), 0.5, 0.25), x = x)
  zk <- zi_cpb(y ~ x, d, zero = ~ 1)
  mk <- underdisp:::.cpb_mean(zk$lambda_full, zk$alpha, FALSE, "ceiling")
  expect_equal(predict(zk, type = "intensity"), mk)
  expect_equal(fitted(zk), (1 - zk$pi_full) * mk)
  zf <- zi_cpb(y ~ x, d, zero = ~ 1, support = "floor")
  expect_equal(predict(zf, type = "intensity"), zf$lambda_full)
})

test_that("the support carries through the hurdle, simulate, and scoring", {
  set.seed(35); n <- 500; x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(y = rhurdle_cpb(n, exp(1.1 + 0.4 * x), 0.5, plogis(0.2 + 0.6 * z)), x = x, z = z)
  hk <- hurdle_cpb(y ~ x, d, participation = ~ z)
  hf <- hurdle_cpb(y ~ x, d, participation = ~ z, support = "floor")
  expect_identical(hk$support, "ceiling"); expect_identical(hf$support, "floor")
  expect_identical(hk$intensity$support, "ceiling")
  ## the in-sample log score is -logLik / n for the intensity's own support
  sk <- score(hk$intensity)
  expect_equal(unname(sk["logscore"]), -hk$intensity$loglik / hk$intensity$n, tolerance = 1e-8)
  sims <- simulate(hk$intensity, nsim = 20, seed = 1)
  top <- ceiling(hk$intensity$ceiling - 1e-9)
  expect_true(all(as.matrix(sims) <= top))
})

test_that("a fit at alpha -> 0 is the two-point distribution of a two-valued series, and print notes the limit", {
  y <- c(rep(4L, 47), rep(5L, 16))                      # a series on two neighbouring counts (Israel's rivalries)
  f <- suppressWarnings(cpb(y ~ 1, data = data.frame(y = y), truncated = TRUE))
  expect_lt(f$alpha, 1e-3)
  expect_lt(abs(f$loglik - (47 * log(47 / 63) + 16 * log(16 / 63))), 1e-3)   # the empirical two-point fit, approached
  expect_output(print(f), "lower limit")
  expect_output(print(summary(f)), "lower limit")
  ci <- underdisp:::.cpb_alpha_ci(f)                     # the interval holds the estimate at the lower limit
  expect_true(ci[["lower"]] <= f$alpha && f$alpha <= ci[["upper"]] && ci[["boundary"]] == 1)
  set.seed(2); fp <- suppressWarnings(cpb(y ~ 1, data = data.frame(y = rpois(200, 3)), truncated = FALSE))
  cp <- underdisp:::.cpb_alpha_ci(fp)                    # and at the Poisson limit, above 0.999
  expect_gt(fp$alpha, 0.999); expect_gte(cp[["upper"]], fp$alpha)
  g <- suppressWarnings(cpb(y ~ 1, data = data.frame(y = c(1L, 2L, 2L, 3L, 1L, 2L, 4L, 2L, 3L, 2L)), truncated = TRUE))
  expect_false(any(grepl("lower limit", capture.output(print(g)))))
})
