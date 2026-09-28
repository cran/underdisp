## dispersion_test() and dispersion_profile()

set.seed(2); n <- 400; x <- rnorm(n)
d_gc <- data.frame(y = rgammacount(n, exp(1 + 0.4 * x), alpha = 2), x = x)
d_po <- data.frame(y = rpois(n, exp(1 + 0.4 * x)), x = x)

test_that("the LR test equals twice the log-likelihood gap to the Poisson refit", {
  m <- count_reg(y ~ x, d_gc, family = "gammacount")
  g <- glm(y ~ x, d_gc, family = poisson())
  dt <- dispersion_test(m)
  expect_identical(dt$calibration, "asymptotic")
  expect_equal(unname(dt$statistic), 2 * (m$loglik - as.numeric(logLik(g))), tolerance = 1e-4)
  expect_equal(dt$p.value, pchisq(unname(dt$statistic), 1, lower.tail = FALSE), tolerance = 1e-8)
  du <- dispersion_test(m, alternative = "under")
  expect_equal(du$p.value, pnorm(sqrt(unname(dt$statistic)), lower.tail = FALSE), tolerance = 1e-8)
  expect_lt(du$p.value, 0.01)                                # regular timing is detected
  do <- dispersion_test(m, alternative = "over")
  expect_gt(do$p.value, 0.5)
  expect_s3_class(dt, "htest")
  expect_output(print(dt), "Likelihood-ratio")
})

test_that("boundary nulls bootstrap by default; B = 0 gives the one-sided Self-Liang mixture", {
  f <- cpb(y ~ x, d_gc, truncated = FALSE, se = "none")
  dt <- dispersion_test(f, B = 0)
  expect_true(dt$boundary)
  expect_false(dt$asymptotic_ok)                       # a boundary null bootstraps by default
  expect_identical(dt$calibration, "asymptotic")
  expect_equal(dt$p.value, 0.5 * pchisq(unname(dt$statistic), 1, lower.tail = FALSE))
  expect_match(paste(dt$note, collapse = " "), "boundary null")
  expect_warning(dispersion_test(f, alternative = "over", B = 0), "one-sided")
  nb <- count_reg(y ~ x, d_po, family = "negbin")
  dn <- dispersion_test(nb, B = 0)
  expect_true(dn$boundary); expect_match(dn$alternative, "over"); expect_false(dn$asymptotic_ok)
  set.seed(7); dnb <- dispersion_test(nb, B = 9)
  expect_identical(dnb$calibration, "parametric bootstrap")
  m <- count_reg(y ~ x, d_gc, family = "gammacount")
  expect_true(dispersion_test(m, B = 0)$asymptotic_ok)  # an interior null keeps the asymptotic distribution
})

test_that("the calibration rule keeps the asymptotic distribution only for small first-order shifts", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  s1 <- underdisp:::.disp_first_order_size
  expect_equal(s1(2, 400, "under"), pnorm(qnorm(0.95) - 2 / sqrt(800), lower.tail = FALSE))
  expect_lt(s1(2, 400, "under"), 0.06)                  # one covariate on 400 rows: asymptotic
  expect_gt(s1(6, 200, "under"), 0.06)                  # five covariates on 200 rows: bootstrap
  expect_lt(s1(6, 200, "over"), 0.05)                   # the shift makes the over test conservative
  set.seed(5); xx <- matrix(rnorm(100 * 5), 100, 5)
  dm <- data.frame(y = rgammacount(100, exp(1 + 0.3 * xx[, 1]), alpha = 1.5), xx)
  m <- count_reg(y ~ ., dm, family = "gammacount")
  set.seed(11); db <- dispersion_test(m, B = 9)
  expect_identical(db$calibration, "parametric bootstrap")
  expect_false(db$asymptotic_ok)
  expect_gte(db$B_ok, 8L)
  expect_true(db$p.value >= 1 / 10 && db$p.value <= 1)
  set.seed(11); expect_identical(dispersion_test(m, B = 9)$p.value, db$p.value)
  da <- dispersion_test(m, B = 0)
  expect_identical(da$calibration, "asymptotic")
  expect_match(paste(da$note, collapse = " "), "first-order size")
  expect_error(dispersion_test(m, B = 2.5), "whole number")
})

test_that("the bootstrap refits reproduce each fit on its own response", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  f <- cpb(y ~ x, d_gc, truncated = FALSE, se = "none")
  expect_equal(underdisp:::.disp_alt_fit(f, f$Y)[["loglik"]], f$loglik, tolerance = 1e-6)
  g <- gec(y ~ x, d_gc, se = "none")
  expect_equal(underdisp:::.disp_alt_fit(g, g$Y)[["loglik"]], g$loglik, tolerance = 1e-6)
  m <- count_reg(y ~ x, d_gc, family = "gammacount")
  expect_equal(underdisp:::.disp_alt_fit(m, m$Y)[["loglik"]], m$loglik, tolerance = 1e-6)
  gl <- glm(y ~ x, d_gc, family = poisson())
  n0 <- underdisp:::.disp_null_fit(f)
  expect_equal(n0$loglik, as.numeric(logLik(gl)), tolerance = 1e-6)
  expect_equal(n0$rate, unname(fitted(gl)), tolerance = 1e-4)
  zt <- d_gc[d_gc$y > 0, ]
  ft <- cpb(y ~ x, zt, truncated = TRUE, se = "none")
  expect_equal(underdisp:::.disp_alt_fit(ft, ft$Y)[["loglik"]], ft$loglik, tolerance = 1e-6)
  set.seed(6); ys <- underdisp:::.disp_simulate(rep(0.3, 2000), TRUE)
  expect_true(all(ys >= 1))
  expect_equal(mean(ys), 0.3 / (1 - exp(-0.3)), tolerance = 0.05)
})

test_that("the fixed-effects null carries the same unit effects and needs the bootstrap", {
  skip_on_cran()
  set.seed(4)
  d <- do.call(rbind, lapply(1:20, function(i) { xx <- rnorm(12); data.frame(unit = i, x = xx,
    y = rcpb(12, exp(rnorm(1, 0.5, 0.4) + 0.4 * xx), 0.6)) }))
  fe <- cpb_fe(y ~ x, d, fe = "unit")
  dt <- dispersion_test(fe, B = 0)
  g <- glm(y ~ x + factor(unit), d, family = poisson())
  expect_equal(unname(dt$loglik["poisson"]), as.numeric(logLik(g)), tolerance = 1e-4)
  expect_true(is.na(dt$p.value))
  expect_match(paste(dt$note, collapse = " "), "fixed effects")
  expect_equal(underdisp:::.disp_null_fit(fe)$rate, unname(fitted(g)), tolerance = 1e-3)
  expect_equal(underdisp:::.disp_alt_fit(fe, fe$Y)[["loglik"]], fe$loglik, tolerance = 1e-6)
  set.seed(12); db <- dispersion_test(fe, B = 19)
  expect_identical(db$calibration, "parametric bootstrap")
  expect_lte(db$p.value, 0.05)
  ge <- gec_fe(y ~ x, d, fe = "unit")
  expect_equal(unname(dispersion_test(ge, B = 0)$loglik["poisson"]), as.numeric(logLik(g)), tolerance = 1e-4)
})

test_that("the auxiliary regression test reproduces the by-hand Cameron-Trivedi statistic", {
  p <- count_reg(y ~ x, d_gc, family = "poisson")
  dt <- dispersion_test(p, method = "auxiliary", alternative = "under")
  mu <- fitted(p); zz <- ((d_gc$y - mu)^2 - d_gc$y) / mu
  t_hand <- summary(lm(zz ~ mu + 0))$coefficients[1, 3]
  expect_equal(unname(dt$statistic), t_hand, tolerance = 1e-10)
  expect_equal(dt$p.value, pnorm(t_hand))
  expect_error(dispersion_test(p), "no dispersion parameter")
})

test_that("the dispersion profile recovers a constant ratio for CPB data", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  set.seed(9); m <- 800; xx <- rnorm(m)
  d <- data.frame(y = rcpb(m, exp(1.5 + 0.5 * xx), 0.5), x = xx)
  f <- cpb(y ~ x, d, truncated = FALSE, se = "none")
  pr <- dispersion_profile(f, negbin = count_reg(y ~ x, d, family = "negbin"), plot = FALSE)
  expect_s3_class(pr, "dispersion_profile")
  expect_equal(nrow(pr), 10L)
  expect_true(all(c("ratio_empirical", "ratio_cpb", "ratio_negbin") %in% names(pr)))
  expect_equal(median(pr$ratio_empirical), 0.5, tolerance = 0.25)
  expect_true(all(pr$ratio_cpb < 0.75))
  expect_output(print(pr), "Dispersion profile")
  pdf(NULL); on.exit(dev.off())
  expect_invisible(plot(pr))
})

test_that("a CPB stopped by the support guard reports the boundary value; a shortfall has no p-value", {
  set.seed(3); xx <- rnorm(300)
  pd <- data.frame(y = rpois(300, exp(1 + 0.3 * xx)), x = xx)
  f <- suppressWarnings(cpb(y ~ x, pd, truncated = FALSE, se = "none", max.support = 40))
  expect_true(f$support_binding)
  expect_lt(f$loglik, f$loglik.null)                    # the guard holds the fit below its Poisson limit
  dt <- dispersion_test(f, B = 0)
  expect_equal(unname(dt$statistic), 0)
  expect_equal(dt$p.value, 0.5)
  expect_match(paste(dt$note, collapse = " "), "boundary value 0")
  expect_output(print(dt), "boundary value 0")
  so <- summary(f)
  expect_equal(so$LR, 0); expect_equal(so$LR.p, 0.5)
  expect_output(print(so), "boundary value 0")
  sh <- underdisp:::.disp_negative_lr(-2, binding = FALSE, "GEC")
  expect_true(sh$shortfall); expect_match(sh$note, "fell short")
})
