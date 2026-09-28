## Regression tests for the round-1 referee findings (statistics and numerics).

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

set.seed(31); n <- 400; x <- rnorm(n); z <- rnorm(n)
d0 <- data.frame(y = rcpb(n, exp(1.2 + 0.5 * x), 0.5), x = x, z = z, u = factor(sample(1:15, n, TRUE)))

test_that("the CPB reports the exact mean of its fitted distribution, not the rate", {
  ## a binomial with non-integer N renormalized over 0..floor(N) has mean != lambda (WSK 1995)
  f <- cpb(y ~ x, d0, truncated = FALSE, se = "none")
  hand <- vapply(f$rate, function(l) { p <- underdisp:::.cpb_pmf_core(l, f$alpha); sum((seq_along(p) - 1) * p) }, numeric(1))
  expect_equal(fitted(f), hand, tolerance = 1e-10)
  expect_equal(mean(fitted(f)), mean(d0$y), tolerance = 0.03)
  expect_equal(predict(f, newdata = d0, type = "rate"), f$rate, tolerance = 1e-10)
  expect_true(any(abs(fitted(f) - f$rate) > 1e-6))
  ## zero-truncated: E(Y | Y >= 1)
  dp <- d0[d0$y > 0, ]; ft <- cpb(y ~ x, dp, truncated = TRUE, se = "none")
  expect_equal(mean(fitted(ft)), mean(dp$y), tolerance = 0.03)
  fd <- first_difference(f, "x", -1, 1)
  expect_equal(fd$from, underdisp:::.cpb_mean(exp(sum(c(1, -1) * coef(f))), f$alpha, FALSE), tolerance = 1e-8)
  ## the hurdle intensity is E(Y | Y > 0)
  h <- hurdle_cpb(y ~ x, d0, participation = ~ z, se = "none")
  expect_equal(fitted(h), h$p_full * underdisp:::.cpb_ztmean(h$lambda_full, h$intensity$alpha), tolerance = 1e-10)
  expect_equal(unname(predict(h, newdata = d0)), unname(fitted(h)), tolerance = 1e-8)
  ## zero-truncated count_reg reports E(Y | Y >= 1)
  ct <- count_reg(y ~ x, dp, family = "poisson", truncated = TRUE)
  mu <- exp(ct$linear.predictors)
  expect_equal(fitted(ct), mu / (1 - exp(-mu)), tolerance = 1e-10)
  expect_equal(mean(fitted(ct)), mean(dp$y), tolerance = 0.03)
})

test_that("the alpha profile interval is re-maximized at each alpha (lower limit below the estimate)", {
  set.seed(32); m <- 300; xx <- rnorm(m)
  d <- data.frame(y = rcpb(m, exp(1 + 0.3 * xx), 0.45), x = xx)
  f <- cpb(y ~ x, d, truncated = FALSE, se = "none")
  ci <- alpha_confint(f)
  expect_lt(ci[[1]], f$alpha - 0.005)
  expect_gt(ci[[2]], f$alpha + 0.005)
  expect_false(attr(ci, "boundary"))
  expect_true(ci[[1]] < 0.45 && 0.45 < ci[[2]])
})

test_that("score() evaluates held-out counts above the in-sample maximum at their own value", {
  d <- data.frame(y = rpois(400, 3))
  f <- count_reg(y ~ 1, d, family = "poisson")
  mu <- fitted(f)[1]
  for (yy in c(14, 29, 200))
    expect_equal(score(f, newdata = data.frame(y = yy))[["logscore"]], -log(max(dpois(yy, mu), 1e-12)), tolerance = 1e-6)
})

test_that("the concentrated fixed-effects likelihood attains the full-dummy maximum", {
  set.seed(33)
  d <- do.call(rbind, lapply(1:12, function(i) { xx <- rnorm(15); data.frame(unit = i, x = xx,
    y = rcpb(15, exp(rnorm(1, 0.5, 0.4) + 0.4 * xx), 0.5)) }))
  fe <- cpb_fe(y ~ x, d, fe = "unit")
  dm <- cpb(y ~ x + factor(unit), d, truncated = FALSE, se = "none")
  expect_gte(fe$loglik, dm$loglik - 0.02)
  ## brute-force check of one unit's intercept at the fitted (beta, alpha)
  u <- 2; rows <- d$unit == u; b <- coef(fe)[["x"]]; a <- fe$alpha
  ll_u <- function(int) sum(log(underdisp:::.cpb_prob_at(exp(int + b * d$x[rows]), a, d$y[rows], FALSE, fe$max.support)))
  grid <- seq(fe$fe[u] - 2, fe$fe[u] + 2, by = 0.002)
  expect_gte(ll_u(fe$fe[u]), max(vapply(grid, ll_u, numeric(1))) - 0.02)
})

test_that("max.support that binds is reported and a guard-stopped LR takes its boundary value", {
  set.seed(801); m <- 300; xx <- rnorm(m)
  big <- data.frame(y = rpois(m, exp(4.6 + 0.2 * xx)), x = xx)      # mean ~ 100, equidispersed
  expect_warning(f <- cpb(y ~ x, big, truncated = FALSE, se = "none", max.support = 500), "max.support")
  expect_true(f$support_binding)
  expect_error(cpb(y ~ x, big, truncated = FALSE, se = "none", max.support = 100), "no feasible fit")
  over <- data.frame(y = rnbinom(m, mu = exp(1.2 + 0.5 * xx), size = 1), x = xx)
  so <- summary(suppressWarnings(cpb(y ~ x, over, truncated = FALSE, se = "none")))
  if (isTRUE(so$support_binding)) { expect_gte(so$LR, 0); expect_lte(so$LR.p, 0.5) } else if (so$LR < 0) expect_true(is.na(so$LR.p)) else expect_true(is.finite(so$LR.p))
  expect_output(print(so), "Poisson")
})

test_that("integer ceilings are feasible in the C++ and R code paths alike", {
  X <- matrix(1, 1, 1); off <- 0
  b <- log(5)                                                   # exp(log(5))/0.5 = 9.999999999999998
  v <- underdisp:::cpb_nll_cpp(c(b, qlogis(0.5)), X, 10L, off, 500L, FALSE)
  expect_true(v < 1e9)
  expect_equal(exp(-v), dcpb(10, 5, 0.5), tolerance = 1e-8)
})

test_that("the pmf consumers carry the fit's offset: in-sample log score equals -logLik/n", {
  E <- runif(n, 0.5, 3); dE <- data.frame(y = rcpb(n, exp(0.9 + 0.4 * x) * E, 0.5), x = x, z = z, logE = log(E))
  f <- cpb(y ~ x, dE, truncated = FALSE, se = "none", offset = "logE")
  expect_equal(score(f)[["logscore"]], -f$loglik / n, tolerance = 1e-8)
  g <- gec(y ~ x, dE, se = "none", offset = "logE")
  expect_equal(score(g)[["logscore"]], -g$loglik / n, tolerance = 1e-8)
  h <- hurdle_count(y ~ x, dE, family = "poisson", participation = ~ z, offset = "logE")
  expect_equal(score(h)[["logscore"]], -h$loglik / n, tolerance = 1e-6)
  hc <- hurdle_cpb(y ~ x, dE, participation = ~ z, offset = "logE", se = "none")
  expect_equal(score(hc)[["logscore"]], -as.numeric(logLik(hc)) / n, tolerance = 1e-6)
})

test_that("the Katz recursion reaches far-tail counts of an unbounded member", {
  expect_equal(dgec(20, 1, 1), dpois(20, 1), tolerance = 1e-10)
  expect_equal(dgec(25, 2, 1.5), dnbinom(25, mu = 2, size = 2 / 0.5), tolerance = 1e-10)
  X <- matrix(1, 5, 1); Y <- c(0L, 1L, 1L, 2L, 30L)
  expect_lt(underdisp:::gec_nll_cpp(c(log(1.5), 0), X, Y, rep(0, 5), 500L), 1e9)
})

test_that("column scaling and restarts make the fits invariant to units and contrasts", {
  d1 <- d0; d1$x <- d1$x * 1000
  for (fam in c("compois", "gammacount", "mpcmp")) {
    a <- count_reg(y ~ x, d0, family = fam); b <- count_reg(y ~ x, d1, family = fam)
    expect_equal(a$loglik, b$loglik, tolerance = 1e-5)
    expect_equal(unname(coef(b)[2]) * 1000, unname(coef(a)[2]), tolerance = 1e-3)
    expect_equal(sqrt(diag(vcov(b)))[2] * 1000, sqrt(diag(vcov(a)))[2], tolerance = 1e-2)
  }
  a <- cpb(y ~ x, d0, truncated = FALSE, se = "none"); b <- cpb(y ~ x, d1, truncated = FALSE, se = "none")
  expect_equal(a$loglik, b$loglik, tolerance = 1e-5)
  g <- factor(sample(c("a", "b", "c"), n, TRUE))
  d2 <- data.frame(y = rcpb(n, exp(0.8 + 0.35 * x + 0.4 * (g == "b") - 0.3 * (g == "c")), 0.5), x = x, g = g)
  fitwith <- function(ctr) { op <- options(contrasts = c(ctr, "contr.poly")); on.exit(options(op))
    cpb(y ~ x + g, d2, truncated = FALSE, se = "none")$loglik }
  ## the pooled CPB likelihood is discontinuous (see ?cpb): two codings can settle on
  ## different teeth, with log-likelihoods differing by the order of the jumps
  expect_lt(abs(fitwith("contr.treatment") - fitwith("contr.sum")), 1.5)
})

test_that("cluster bootstraps refuse a single cluster and first differences refuse interactions", {
  skip_on_cran()
  d1 <- transform(d0, one = 1L)
  expect_error(cpb(y ~ x, d1, truncated = FALSE, se = "bootstrap", B = 5, cluster = "one"), "at least two clusters")
  expect_error(gec(y ~ x, d1, se = "bootstrap", B = 5, cluster = "one"), "at least two clusters")
  fi <- cpb(y ~ x * z, d0, truncated = FALSE, se = "none")
  expect_error(first_difference(fi, "x", 0, 1), "interaction")
  ci <- count_reg(y ~ x * z, d0, family = "poisson")
  expect_error(first_difference(ci, "x", 0, 1), "interaction")
})

test_that("first differences hold an offset at its mean", {
  E <- runif(n, 0.5, 3); dE <- data.frame(y = rcpb(n, exp(0.9 + 0.4 * x) * E, 0.5), x = x, logE = log(E))
  f <- cpb(y ~ x, dE, truncated = FALSE, se = "none", offset = "logE")
  fd <- first_difference(f, "x", 0, 0)
  expect_equal(fd$from, underdisp:::.cpb_mean(exp(mean(log(E)) + coef(f)[[1]]), f$alpha, FALSE), tolerance = 1e-8)
})
