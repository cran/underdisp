## Split-panel jackknife on the concentrated FE estimators, INCLUDING the
## validity gate. The assumption-violating scenarios are tested explicitly:
## a time-homogeneous panel must be corrected, a dispersion-regime-change panel
## must be REFUSED (the Dhaene-Jochmans time-homogeneity requirement), and
## short panels must warn. DGP mirrors the committed fixed-effects bias Monte
## Carlo (main_fe_bias_mc.R).

make_panel <- function(seed, N = 40, TT = 8, alpha = 0.5, beta = 0.3) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    x <- rnorm(TT); lam <- exp(1.0 + rnorm(1, 0, 0.5) + beta * x)
    data.frame(unit = i, tt = seq_len(TT), x = x, y = rcpb(TT, lam, alpha))
  }))
}

## alpha shifts .3 -> .7 at mid-panel: the two halves estimate different
## parameters, so the split-panel identity fails and the gate must refuse
make_regime_panel <- function(seed, N = 40, TT = 12) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    x <- rnorm(TT); lam <- exp(1.0 + rnorm(1, 0, 0.4) + 0.3 * x)
    a <- ifelse(seq_len(TT) <= TT / 2, 0.3, 0.7)
    y <- vapply(seq_len(TT), function(t) rcpb(1, lam[t], a[t]), integer(1))
    data.frame(unit = i, tt = seq_len(TT), x = x, y = y)
  }))
}

test_that("cpb_fe jackknife corrects on a time-homogeneous panel (audit trail intact)", {
  skip_on_cran()
  d <- make_panel(101)
  raw <- cpb_fe(y ~ x, data = d, fe = "unit")
  jk  <- cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt")
  expect_identical(jk$bias_correct, "jackknife")        # gate did NOT refuse
  expect_identical(jk$uncorrected$alpha, raw$alpha)
  expect_equal(jk$uncorrected$coefficients, raw$coefficients)
  expect_lt(abs(jk$alpha - 0.5), abs(raw$alpha - 0.5))  # moves toward truth
  expect_lt(abs(jk$alpha - 0.5), 0.12)
  expect_true(all(is.finite(jk$fitted.values)))
  expect_output(print(jk), "jackknife bias-corrected")
  s <- simulate(jk, nsim = 3, seed = 2)                 # DHARMa path intact
  expect_identical(dim(as.matrix(s)), c(nrow(d), 3L))
})

test_that("gate REFUSES on a dispersion regime change (time-heterogeneity)", {
  skip_on_cran()
  d <- make_regime_panel(202)
  expect_warning(jk <- cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt"),
                 "REFUSED")
  expect_identical(jk$bias_correct, "none")             # ML fit returned
  expect_null(jk$uncorrected)
  raw <- cpb_fe(y ~ x, data = d, fe = "unit")
  expect_equal(jk$alpha, raw$alpha, tolerance = 1e-8)   # identical to uncorrected
})

test_that("gate REFUSES on a smooth unmodeled trend (residual time-trend check)", {
  skip_on_cran()
  set.seed(203)
  d <- do.call(rbind, lapply(seq_len(40), function(i) {
    TT <- 10L; x <- rnorm(TT)
    lam <- exp(1.0 + rnorm(1, 0, 0.4) + 0.3 * x + log(3) * (seq_len(TT) - 1) / (TT - 1))
    data.frame(unit = i, tt = seq_len(TT), x = x, y = rcpb(TT, lam, 0.5))
  }))
  expect_warning(jk <- cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt"),
                 "residuals trend with time")
  expect_identical(jk$bias_correct, "none")
})

test_that("gec_fe jackknife corrects delta and shares the gate", {
  skip_on_cran()
  d <- make_panel(303)
  raw <- gec_fe(y ~ x, data = d, fe = "unit")
  jk  <- gec_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt")
  expect_identical(jk$bias_correct, "jackknife")
  expect_identical(jk$uncorrected$delta, raw$delta)
  expect_lte(abs(jk$delta - 0.5), abs(raw$delta - 0.5) + 0.02)
  expect_output(print(jk), "jackknife bias-corrected")
  d2 <- make_regime_panel(304)
  expect_warning(g2 <- gec_fe(y ~ x, data = d2, fe = "unit", bias_correct = "jackknife", time = "tt"),
                 "REFUSED")
  expect_identical(g2$bias_correct, "none")
})

test_that("hurdle_cpb passes bias_correct through to the FE intensity", {
  skip_on_cran()
  d <- make_panel(405, N = 30, TT = 14)
  d$y[sample(nrow(d), nrow(d) %/% 5)] <- 0L
  fit <- hurdle_cpb(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt")
  expect_identical(fit$intensity$bias_correct, "jackknife")
  expect_false(is.null(fit$intensity$uncorrected))
})

test_that("short panels warn; zi_cpb no longer offers bias_correct", {
  skip_on_cran()
  d <- make_panel(506, N = 30, TT = 4)
  expect_warning(cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt"),
                 "median panel length")
  dz <- make_panel(507, N = 20, TT = 8)
  dz$y[sample(nrow(dz), 40)] <- 0L
  expect_error(zi_cpb(y ~ x, dz, zero = ~ 1, fe = "unit", bias_correct = "jackknife"),
               "unused argument")
})

test_that("the jackknife needs `time`, and splits each unit's series in time order whatever the row order", {
  skip_on_cran()
  d <- make_regime_panel(202)
  expect_error(cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife"), "needs 'time'")
  expect_error(gec_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife"), "needs 'time'")
  set.seed(9); sh <- d[sample(nrow(d)), ]
  ## the regime change is refused however the rows are supplied
  expect_warning(jk <- cpb_fe(y ~ x, data = sh, fe = "unit", bias_correct = "jackknife", time = "tt"), "REFUSED")
  expect_identical(jk$bias_correct, "none")
  d2 <- make_panel(101); set.seed(10); sh2 <- d2[sample(nrow(d2)), ]
  a1 <- cpb_fe(y ~ x, data = d2, fe = "unit", bias_correct = "jackknife", time = "tt")
  a2 <- cpb_fe(y ~ x, data = sh2, fe = "unit", bias_correct = "jackknife", time = "tt")
  expect_equal(a1$alpha, a2$alpha, tolerance = 1e-8)
  expect_equal(unname(fitted(a2)), unname(fitted(a1))[match(rownames(sh2), rownames(d2))], tolerance = 1e-8)
  expect_error(cpb_fe(y ~ x, data = d2, fe = "unit", time = "nope"), "'time' must name a column")
})

test_that("a jackknife fit's bootstrap corrects every replicate on King's support, not on the floor's", {
  skip_on_cran()
  d <- make_panel(101, N = 30, TT = 8)
  set.seed(5); f0 <- suppressMessages(cpb_fe(y ~ x, data = d, fe = "unit", se = "bootstrap", B = 12))
  set.seed(5); f1 <- suppressMessages(cpb_fe(y ~ x, data = d, fe = "unit", se = "bootstrap", B = 12,
                                             bias_correct = "jackknife", time = "tt"))
  expect_identical(f1$bias_correct, "jackknife")
  expect_null(getOption("underdisp.jackknife_replicate"))            # the replicate mode does not leak
  ok <- stats::complete.cases(f0$boot) & stats::complete.cases(f1$boot)
  expect_gt(sum(ok), 8)
  expect_gt(mean(f1$boot[ok, "alpha"]), mean(f0$boot[ok, "alpha"]))  # the same resamples, corrected
  expect_false(isTRUE(all.equal(f1$se.beta, f0$se.beta)))
  expect_equal(unname(f1$ci.beta[, "lower"]), unname(f1$coefficients - stats::qnorm(0.975) * f1$se.beta))
  ## the floor support keeps uncorrected replicates, the interval recentred on the corrected estimate
  set.seed(7); g0 <- suppressWarnings(suppressMessages(cpb_fe(y ~ x, data = d, fe = "unit", se = "bootstrap", B = 8,
                                                              support = "floor")))
  set.seed(7); g1 <- suppressWarnings(suppressMessages(cpb_fe(y ~ x, data = d, fe = "unit", se = "bootstrap", B = 8,
                                                              support = "floor", bias_correct = "jackknife", time = "tt")))
  expect_equal(g1$boot, g0$boot)
  expect_equal(g1$se.beta, g0$se.beta)
})
