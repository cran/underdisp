## Regression tests for the round-2 referee items (see
## verification/referee_panels_2026-09/round2/response_r2.md): each block guards
## one item by the behaviour the referee reproduced.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

set.seed(21); n <- 240; x <- rnorm(n); z <- rnorm(n); u <- rep(1:20, each = 12)
d  <- data.frame(y = rcpb(n, exp(1.2 + 0.4 * x + rnorm(20, 0, 0.3)[u]), 0.5), x = x, z = z, unit = u)
dz <- d; dz$y <- ifelse(runif(n) < plogis(-1 + 0.8 * z), 0L, dz$y)
dp <- d[d$y > 0, ]

test_that("zero-truncated fits are scored with their own predictive distribution (A-R1)", {
  g <- gec(y ~ x, dp, truncated = TRUE, se = "none")
  expect_equal(unname(score(g)["logscore"]), -as.numeric(logLik(g)) / nobs(g), tolerance = 1e-10)
  for (f in c("poisson", "negbin", "gammacount")) {
    cr <- count_reg(y ~ x, dp, family = f, truncated = TRUE)
    expect_equal(unname(score(cr)["logscore"]), -as.numeric(logLik(cr)) / nobs(cr), tolerance = 1e-10, label = f)
    expect_equal(unname(score(cr, newdata = dp)["logscore"]), unname(score(cr)["logscore"]), tolerance = 1e-10)
  }
  expect_equal(unname(score(g, newdata = dp)["logscore"]), unname(score(g)["logscore"]), tolerance = 1e-10)
})

test_that("a misaligned or incomplete vector cluster is refused by every estimator (A-R2, A-R12)", {
  bad3 <- c("a", "b", "c"); clna <- as.character(d$unit); clna[1:5] <- NA
  expect_error(cpb_fe(y ~ x, d, fe = "unit", cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(gec_fe(y ~ x, d, fe = "unit", cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(hurdle_cpb(y ~ x, dz, participation = ~ z, cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(hurdle_gec(y ~ x, dz, participation = ~ z, cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(zi_cpb(y ~ x, dz, zero = ~ z, cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(zi_gec(y ~ x, dz, zero = ~ z, cluster = bad3, se = "bootstrap", B = 3), "one value per row")
  expect_error(cpb_fe(y ~ x, d, fe = "unit", cluster = clna, se = "bootstrap", B = 3), "missing values")
  expect_error(gec_fe(y ~ x, d, fe = "unit", cluster = clna, se = "bootstrap", B = 3), "missing values")
  ## a correctly aligned vector works where a column name does, in hurdle_count too
  clv <- as.character(dz$unit)
  h1 <- hurdle_count(y ~ x, dz, family = "poisson", participation = ~ z, cluster = clv, se = "cluster")
  h2 <- hurdle_count(y ~ x, dz, family = "poisson", participation = ~ z, cluster = "unit", se = "cluster")
  expect_equal(h1$intensity$se.beta, h2$intensity$se.beta)
  ## supplying cluster selects the cluster-robust SE in hurdle_count as in its siblings (A-R7)
  h3 <- hurdle_count(y ~ x, dz, family = "poisson", participation = ~ z, cluster = "unit")
  expect_identical(h3$intensity$se.type, "cluster")
})

test_that("a fractional response is refused by the zero-inflated mixtures (A-R3)", {
  df <- dz; df$y <- df$y + 0.4
  expect_error(zi_cpb(y ~ x, df, zero = ~ z, se = "none"), "integer")
  expect_error(zi_gec(y ~ x, df, zero = ~ z, se = "none"), "integer")
})

test_that("fitted(), predict() and the stored fields agree within each class (A-R4, A-R11)", {
  h <- hurdle_gec(y ~ x, dz, participation = ~ z, se = "none")
  expect_equal(unname(fitted(h)), unname(predict(h, type = "response")), tolerance = 1e-12)
  cr <- count_reg(y ~ x, dp, family = "gammacount", truncated = TRUE)
  expect_equal(cr$fitted.values, fitted(cr)); expect_equal(cr$residuals, residuals(cr))
  expect_equal(unname(fitted(cr)), unname(predict(cr, type = "response")), tolerance = 1e-12)
  expect_true(all(cr$mu > 0) && !isTRUE(all.equal(cr$mu, cr$fitted.values)))   # natural parameter kept apart
  zc <- zi_count(y ~ x, dz, family = "poisson", zero = ~ z)
  expect_equal(zc$fitted.values, fitted(zc))
  expect_equal(unname(fitted(zc)), unname(predict(zc, type = "response")), tolerance = 1e-12)
})

test_that("summary.zi_count prints the stored inflation SEs (A-R5)", {
  zc <- zi_count(y ~ x, dz, family = "poisson", zero = ~ z)
  s <- capture.output(print(summary(zc)))
  expect_true(sum(grepl("Std. Error", s)) >= 2)
})

test_that("two-part coefficient names follow one rule across the methods (A-R6)", {
  fits <- list(hurdle_cpb(y ~ x, dz, participation = ~ z, se = "none"),
               hurdle_gec(y ~ x, dz, participation = ~ z, se = "none"),
               hurdle_count(y ~ x, dz, family = "poisson", participation = ~ z),
               zi_count(y ~ x, dz, family = "poisson", zero = ~ z))
  skip_if_not_installed("broom")
  for (f in fits) {
    td <- broom::tidy(f)
    if (!is.null(vcov(f))) expect_setequal(td$term, rownames(vcov(f)))
    if (!inherits(f, c("hurdle_cpb", "hurdle_gec"))) expect_setequal(td$term, rownames(confint(f)))
    expect_true(all(grepl("^(participation|intensity|count|zero):", td$term)))
  }
  skip_if_not_installed("texreg")
  tr <- texreg::extract(fits[[3]])
  expect_setequal(tr@coef.names, broom::tidy(fits[[3]])$term)
})

test_that("the Nelder-Mead polish chains restarts until the gain drops below tol (A-R8)", {
  rosen <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2 + 100 * (p[3] - p[2]^2)^2 + (1 - p[2])^2
  f0 <- stats::optim(c(-1.2, 1, 1), rosen, method = "Nelder-Mead", control = list(maxit = 60))
  p1 <- underdisp:::.ud_nm_polish(f0, rosen, maxit = 60, reltol = 1e-8, restarts = 1L)
  p6 <- underdisp:::.ud_nm_polish(f0, rosen, maxit = 60, reltol = 1e-8, restarts = 6L)
  expect_lt(p6$value, p1$value - 1e-3)
})

test_that("zi_cpb is invariant to the units of a covariate (A-R9)", {
  f1 <- zi_cpb(y ~ x, dz, zero = ~ z, se = "none")
  d2 <- dz; d2$x <- 1000 * d2$x; d2$z <- 100 * d2$z
  f2 <- zi_cpb(y ~ x, d2, zero = ~ z, se = "none")
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-8)
  expect_equal(unname(coef(f1)$count["x"]), 1000 * unname(coef(f2)$count["x"]), tolerance = 1e-6)
  expect_equal(unname(f1$zero_coef["z"]), 100 * unname(f2$zero_coef["z"]), tolerance = 1e-6)
})

test_that("every bootstrapping summary states its survivor count (C-N11)", {
  skip_on_cran()
  for (f in list(cpb(y ~ x, d, truncated = FALSE, se = "bootstrap", B = 4),
                 gec(y ~ x, d, se = "bootstrap", B = 4),
                 cpb_fe(y ~ x, d, fe = "unit", se = "bootstrap", B = 4),
                 gec_fe(y ~ x, d, fe = "unit", se = "bootstrap", B = 4)))
    expect_true(any(grepl("bootstrap resamples converged", capture.output(print(summary(f))))))
})

test_that("rows dropped for a missing unit identifier take vector arguments with them", {
  dn <- d; dn$unit[1:6] <- NA
  off <- rep(0.1, n); wt <- rep(1, n)
  f <- cpb_fe(y ~ x, dn, fe = "unit", offset = off, weights = wt, se = "none")
  expect_equal(f$n, n - 6)
  g <- count_reg(y ~ x, dn, family = "poisson", fe = "unit", offset = off, weights = wt)
  expect_equal(g$n, n - 6)
})

test_that("zi_gec does not depend on the units of a covariate", {
  f1 <- zi_gec(y ~ x, dz, zero = ~ z, se = "none")
  d2 <- dz; d2$x <- 1000 * d2$x; d2$z <- 100 * d2$z
  f2 <- zi_gec(y ~ x, d2, zero = ~ z, se = "none")
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-6)
  expect_equal(unname(coef(f1)$count["x"]), 1000 * unname(coef(f2)$count["x"]), tolerance = 1e-4)
})
