## New data, and the fixed-effects units a fit has not seen or holds at a limit (?score, ?compare_models, Degrees of
## freedom): every two-part class predicts and scores such rows, an unseen unit at the reference unit of the units the
## data identify and a unit at a limit at that limit. Also other paths through the same machinery: several
## fixed-effects columns, the floor support's zi_cpb(fe) df, the ML mixture's reference unit, the support of ame()'s
## refits, the dispersion test's fixed-effects null, the 1e5-term cap, and compare_models()'s note.

## 14 units x 6 periods: u01 (the dummies' reference level) and u07 have no positive count, u02 has every count positive
nd_panel <- function(seed = 1) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(14), function(u) {
    x <- rnorm(6); z <- rnorm(6); a <- rnorm(1, 0.8, 0.3); g <- rnorm(1, 0, 0.5)
    part <- rbinom(6, 1, plogis(g + 1.2 * z))
    y <- ifelse(part == 1, pmax(1L, rbinom(6, pmax(1, round(exp(a + 0.4 * x) / 0.5)), 0.5)), 0L)
    data.frame(unit = sprintf("u%02d", u), t = seq_len(6), x = x, z = z, y = as.integer(y))
  }))
  d$y[d$unit %in% c("u01", "u07")] <- 0L
  d$y[d$unit == "u02" & d$y == 0L] <- 1L
  d
}
q <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("a unit the fit has not seen, or newdata without the unit column, enters at the reference unit", {
  skip_on_cran()
  d <- nd_panel(); nd <- data.frame(unit = "u99", t = 7L, x = 0.3, z = 0, y = 1L)
  fits <- list(
    hurdle_cpb = q(hurdle_cpb(y ~ x, d, participation = ~ z, fe = "unit", part_fe = "unit")),
    hurdle_gec = q(hurdle_gec(y ~ x, d, participation = ~ z, part_fe = "unit")),
    hurdle_count = q(hurdle_count(y ~ x, d, family = "poisson", participation = ~ z, fe = "unit", part_fe = "unit",
                                  se = "none")),
    zi_cpb = q(zi_cpb(y ~ x, d, zero = ~ z, zero_fe = "unit")),
    zi_gec = q(zi_gec(y ~ x, d, zero = ~ z, zero_fe = "unit")),
    zi_count = q(zi_count(y ~ x, d, family = "poisson", zero = ~ z, fe = "unit", zero_fe = "unit", se = "none")))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_true(is.finite(q(predict(f, newdata = nd))), info = nm)
    expect_true(all(is.finite(q(score(f, newdata = nd)))), info = nm)
    expect_equal(q(predict(f, newdata = nd[, c("t", "x", "z", "y")])), q(predict(f, newdata = nd)), info = nm)
  }
  ## the reference is that of the units the data identify (first_difference()'s part_xref), not the mean of every
  ## dummy, which the never- and always-participating units' diverging dummies would pull to either end
  h <- fits$hurdle_cpb; b <- coef(h$participation); dm <- startsWith(names(b), "factor(unit)")
  pref <- plogis(b[["(Intercept)"]] + sum(h$part_xref[names(b)[dm]] * b[dm]))
  expect_equal(predict(h, newdata = nd, type = "participation"), pref)
  expect_equal(1 - underdisp:::.pmf_and_y_newdata(h, nd, kmax = 10L)$P[1, 1], pref)
  ## the floor support scores it as before 0.2.0, at the mean of the factor's estimated effects
  hf <- q(hurdle_cpb(y ~ x, d, participation = ~ z, fe = "unit", part_fe = "unit", support = "floor"))
  bf <- coef(hf$participation); dmf <- startsWith(names(bf), "factor(unit)")
  expect_equal(1 - underdisp:::.pmf_and_y_newdata(hf, nd, kmax = 10L)$P[1, 1],
               plogis(bf[["(Intercept)"]] + mean(c(0, bf[dmf]))))
  ## rows the fit can evaluate are unchanged: a seen unit predicts as the glm does
  ns <- d[d$t == 1, ]
  expect_identical(predict(h, newdata = ns[!ns$unit %in% c("u01", "u07", "u02"), ], type = "participation"),
                   as.numeric(predict(h$participation, newdata = ns[!ns$unit %in% c("u01", "u07", "u02"), ],
                                      type = "response")))
})

test_that("hurdle_count(fe =) predicts and scores the rows of a unit its intensity never saw", {
  d <- nd_panel()
  hc <- function(tr) hurdle_count(y ~ x, tr, family = "poisson", participation = ~ z, fe = "unit", part_fe = "unit",
                                  se = "none")
  cv <- q(cv_score(hc, d, folds = rep(1:3, length.out = nrow(d))))
  expect_false(anyNA(cv$logscore_i))
  f <- q(hc(d)); u7 <- d[d$unit == "u07", ]
  expect_identical(predict(f, newdata = u7), rep(0, nrow(u7)))
  expect_identical(underdisp:::.pmf_and_y_newdata(f, u7, kmax = 5L)$P[, 1], rep(1, nrow(u7)))
})

test_that("the mixtures give a zero_fe unit at a limit its structural-zero probability there", {
  skip_on_cran()
  d <- nd_panel()
  ## covariates far out, where the optimizer's finite dummies would give either unit the opposite probability
  nd <- data.frame(unit = c("u07", "u02"), t = 7L, x = 0, z = c(40, -40), y = 0L)
  for (f in list(q(zi_cpb(y ~ x, d, zero = ~ z, zero_fe = "unit")), q(zi_gec(y ~ x, d, zero = ~ z, zero_fe = "unit")),
                 q(zi_count(y ~ x, d, family = "poisson", zero = ~ z, zero_fe = "unit", se = "none")))) {
    expect_setequal(f$zero_fe_one_units, "u02")
    expect_identical(predict(f, newdata = nd, type = "zero"), c(1, 0))
    expect_identical(underdisp:::.pmf_and_y_newdata(f, nd, kmax = 5L)$P[1, 1], 1)
  }
  ## a count-equation unit with no positive count predicts a mean of zero, its zeros certain
  zc <- q(zi_count(y ~ x, d, family = "poisson", zero = ~ z, fe = "unit", se = "none"))
  expect_identical(predict(zc, newdata = nd[1, ], type = "count"), 0)
  expect_identical(underdisp:::.pmf_and_y_newdata(zc, nd[1, ], kmax = 5L)$P[1, ], c(1, rep(0, 5)))
  zk <- q(zi_cpb(y ~ x, d, zero = ~ z, fe = "unit"))
  expect_identical(predict(zk, newdata = nd[1, ]), 0)
  expect_identical(underdisp:::.pmf_and_y_newdata(zk, nd[1, ], kmax = 5L)$P[1, ], c(1, rep(0, 5)))
  ## the floor support predicts as before 0.2.0
  zf <- q(zi_cpb(y ~ x, d, zero = ~ z, zero_fe = "unit", support = "floor"))
  expect_length(zf$zero_fe_one_units, 0L)
  expect_length(zf$zero_fe_zero_units, 0L)
})

test_that("several fixed-effects columns enter the binary equations, as in 0.1.x", {
  skip_on_cran()
  d <- nd_panel(); d$per <- factor(d$t)
  hc <- q(hurdle_count(y ~ x, d, family = "poisson", participation = ~ z, part_fe = c("unit", "per"), se = "none"))
  expect_setequal(hc$part_zero_units, c("u01", "u07"))
  expect_identical(hc$df, length(coef(hc$participation)) - 2L + hc$intensity$df)
  zc <- q(zi_count(y ~ x, d, family = "poisson", zero = ~ z, zero_fe = c("unit", "per"), se = "none"))
  expect_setequal(zc$zero_fe_zero_units, c("u01", "u07"))
  expect_true(is.finite(zc$loglik))
  for (f in list(q(zi_cpb(y ~ x, d, zero = ~ z, zero_fe = c("unit", "per"))),
                 q(zi_gec(y ~ x, d, zero = ~ z, zero_fe = c("unit", "per"))),
                 q(hurdle_cpb(y ~ x, d, participation = ~ z, part_fe = c("unit", "per"))),
                 q(hurdle_gec(y ~ x, d, participation = ~ z, part_fe = c("unit", "per"))))) {
    expect_true(is.finite(as.numeric(logLik(f))))
    expect_true(is.finite(q(predict(f, newdata = d[1:3, ])))[1])
  }
})

test_that("the floor support's zi_cpb(fe, method = 'ml') counts every unit, as before 0.2.0", {
  d <- nd_panel()
  z <- q(zi_cpb(y ~ x, d, zero = ~ z, fe = "unit", support = "floor"))
  expect_identical(z$df, length(z$coefficients) + length(z$fe_hat) + 1L + length(z$zero_coef))
})

test_that("King's zi_cpb(fe, method = 'ml') averages its reference unit over every unit with a positive count", {
  d <- nd_panel(); one <- d$unit %in% c("u03", "u05"); d$y[one] <- pmin(d$y[one], 1L)
  z <- q(zi_cpb(y ~ x, d, zero = ~ z, fe = "unit"))
  pos <- tapply(d$y > 0, d$unit, any)
  expect_equal(z$int_fe_ref, mean(z$fe_hat[names(pos)[pos]]))
  ze <- q(zi_cpb(y ~ x, d, zero = ~ z, fe = "unit", method = "em"))
  expect_equal(predict(z, newdata = data.frame(unit = "u99", x = 0, z = 0)),
               predict(ze, newdata = data.frame(unit = "u99", x = 0, z = 0)), tolerance = 0.05)
})

test_that("ame(B > 0) refits an object saved before 0.2.0 on the floor support", {
  skip_on_cran()
  set.seed(2); x <- rnorm(150)
  dd <- data.frame(y = rcpb(150, exp(1.2 + 0.3 * x), 0.5, truncated = TRUE, support = "floor"), x = x)
  fo <- q(cpb(y ~ x, dd, support = "floor"))
  old <- fo; old$support <- NULL; old$call$support <- NULL       # as an object saved by 0.1.x
  set.seed(5); a1 <- q(ame(old, "x", dd, from = -1, to = 1, B = 3))
  set.seed(5); a2 <- q(ame(fo, "x", dd, from = -1, to = 1, B = 3))
  expect_equal(as.data.frame(a1), as.data.frame(a2))
})

test_that("dispersion_test()'s fixed-effects null does not depend on centring a covariate", {
  set.seed(4); u <- rep(1:20, each = 8); yr <- rep(2001:2008, 20)
  y <- rpois(160, exp(0.4 + 0.05 * (yr - 2004.5) + rnorm(20, 0, 0.3)[u])); y[u == 1] <- 0L
  d <- data.frame(y = y, u = u, year = yr, yc = yr - mean(yr))
  t1 <- q(dispersion_test(cpb_fe(y ~ year, d, fe = "u"), B = 0))
  t2 <- q(dispersion_test(cpb_fe(y ~ yc, d, fe = "u"), B = 0))
  ## where cpb_fe() stops sets the gap: 4e-7 relative on x86_64, 1.6e-6 on Linux ARM64; 0.1.x was off by a factor of 1e10
  expect_equal(unname(t1$statistic), unname(t2$statistic), tolerance = 1e-4)
})

test_that("past the 1e5-term cap, predicted probabilities and draws are NA with a warning", {
  set.seed(1); x <- rnorm(200)
  d <- data.frame(y = rcpb(200, exp(1 + 0.3 * x), 0.5, truncated = TRUE), x = x)
  f <- q(cpb(y ~ x, d))
  nd <- data.frame(x = (log(2e5) - f$coefficients[1]) / f$coefficients[2])
  expect_warning(p <- predict(f, nd, type = "prob", at = 5), "1e5-term cap")
  expect_true(is.na(p))
  g <- f; g$rate[1] <- 2e5
  expect_warning(s <- simulate(g, nsim = 3, seed = 2), "1e5-term cap")
  expect_true(all(is.na(unlist(s[1, ]))))
  expect_false(anyNA(unlist(s[-1, ])))
})

test_that("compare_models() notes the degrees of freedom a single-equation fixed-effects fit charges", {
  d <- nd_panel()
  expect_message(q2 <- suppressWarnings(compare_models(zi = q(zi_cpb(y ~ x, d, zero = ~ z, fe = "unit")),
                                                       fe = q(cpb_fe(y ~ x, d, fe = "unit")))),
                 "count the effect of each fixed-effects unit with no positive count|counts the effect")
  expect_s3_class(q2, "data.frame")
})
