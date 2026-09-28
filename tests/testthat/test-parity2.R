## Method parity across the eleven model classes (referee A, round 1): every
## class answers confint()/vcov()/irr()/first_difference()/update(), the
## refusals are informative, and the argument checks are uniform.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

set.seed(21); n <- 300; x <- rnorm(n); z <- rnorm(n)
d0 <- data.frame(y = rcpb(n, exp(1.2 + 0.5 * x), 0.5), x = x, z = z, u = factor(sample(1:15, n, TRUE)))

test_that("confint/vcov/irr/first_difference/update work on every class with inference", {
  skip_on_cran()
  fits <- list(cpb = cpb(y ~ x, d0, truncated = FALSE, se = "bootstrap", B = 25),
               cpb_fe = cpb_fe(y ~ x, d0, fe = "u", se = "bootstrap", B = 15),
               gec = gec(y ~ x, d0, se = "bootstrap", B = 25),
               gec_fe = gec_fe(y ~ x, d0, fe = "u", se = "bootstrap", B = 15),
               hurdle_cpb = hurdle_cpb(y ~ x, d0, participation = ~ z, se = "bootstrap", B = 25),
               hurdle_gec = hurdle_gec(y ~ x, d0, participation = ~ z, se = "bootstrap", B = 25),
               zi_cpb = zi_cpb(y ~ x, d0, zero = ~ z, se = "bootstrap", B = 15),
               zi_gec = zi_gec(y ~ x, d0, zero = ~ z, se = "bootstrap", B = 15),
               count_reg = count_reg(y ~ x, d0, family = "negbin"),
               hurdle_count = hurdle_count(y ~ x, d0, family = "poisson", participation = ~ z),
               zi_count = zi_count(y ~ x, d0, family = "poisson", zero = ~ z))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    ci <- confint(f); expect_true(is.matrix(ci) && ncol(ci) == 2L, label = paste(nm, "confint"))
    expect_true(all(ci[, 1] <= ci[, 2], na.rm = TRUE), label = paste(nm, "confint order"))
    v <- vcov(f); expect_true(is.matrix(v), label = paste(nm, "vcov"))
    expect_s3_class(irr(f), "ud_irr")
    expect_s3_class(first_difference(f, "x", 0, 1), "ud_fd")
    f2 <- update(f, . ~ . + z); expect_s3_class(f2, class(f)[1L])
    expect_false(is.null(f$converged), label = paste(nm, "converged stored"))
  }
  expect_equal(nrow(implied_ceiling(fits$cpb_fe)), n)
  expect_equal(nrow(implied_ceiling(hurdle_cpb(y ~ x, d0, participation = ~ z, fe = "u", se = "none"))), n)
  expect_setequal(unique(irr(fits$hurdle_cpb)$equation), c("count", "participation"))
  expect_setequal(unique(irr(fits$zi_cpb)$equation), c("count", "inflation"))
  expect_true(all(is.finite(confint(fits$zi_count))))            # list coef() no longer breaks confint
  expect_equal(nrow(confint(fits$count_reg)), 3L)                # coefficients + dispersion
})

test_that("fits without inference refuse informatively or return point estimates", {
  f <- cpb(y ~ x, d0, truncated = FALSE, se = "none")
  expect_null(vcov(f))
  expect_error(confint(f), "bootstrap")
  expect_equal(unique(irr(f)$method), "none")
  fd <- first_difference(f, "x", 0, 1); expect_equal(fd$method, "none"); expect_true(is.finite(fd$diff))
  g <- gec(y ~ x, d0, se = "none")
  expect_null(vcov(g)); expect_error(confint(g), "No inference")
  fe <- cpb_fe(y ~ x, d0, fe = "u")
  expect_error(confint(fe), "No inference")
  expect_null(vcov(hurdle_cpb(y ~ x, d0, participation = ~ z, se = "none")))
})

test_that("argument checks are uniform across the family", {
  expect_error(count_reg(y ~ x, d0, foo = 1), "unused argument")
  expect_error(hurdle_count(y ~ x, d0, foo = 1), "unused argument")
  expect_error(zi_count(y ~ x, d0, foo = 1), "unused argument")
  expect_warning(cpb(y ~ x, d0, truncated = FALSE, cluster = "u"), "only affects the bootstrap")
  expect_warning(gec(y ~ x, d0, cluster = "u"), "only affects the bootstrap")
  expect_error(cpb_fe(y ~ x, d0, fe = "u", offset = c(0.1, -0.2, 0.05)), "one value per row")
  expect_error(gec_fe(y ~ x, d0, fe = "u", offset = c(0.1, -0.2, 0.05)), "one value per row")
  expect_error(zi_cpb(y ~ x, d0, zero = ~ z, offset = c(0.1, -0.2, 0.05)), "one value per row")
  expect_error(gec(y ~ x, transform(d0, y = y + 0.5), se = "none"), "integer")
  expect_error(gec_fe(y ~ x, transform(d0, y = y + 0.5), fe = "u"), "integer")
  expect_error(cpb_fe(y ~ x, transform(d0, y = y + 0.5), fe = "u", truncated = FALSE), "integer")
  z1 <- zi_cpb(y ~ x, d0, zero = ~ z)
  expect_identical(as.character(z1$call[[1L]]), "zi_cpb")
})

test_that("rank deficiency and within-unit constancy are refused with the column named", {
  expect_error(cpb(y ~ x + I(2 * x), d0, truncated = FALSE, se = "none"), "I\\(2 \\* x\\)")
  expect_error(gec(y ~ x + I(2 * x), d0, se = "none"), "rank deficient")
  expect_error(count_reg(y ~ x + I(2 * x), d0), "rank deficient")
  g <- factor(sample(c("a", "b", "z"), n, TRUE, prob = c(.5, .5, 0)), levels = c("a", "b", "z"))
  ## an unused level is not a rank deficiency: it is dropped from the model frame, as lm() and glm() drop it
  expect_equal(suppressWarnings(cpb(y ~ g, transform(d0, g = g), truncated = FALSE, se = "none"))$loglik,
               suppressWarnings(cpb(y ~ g, transform(d0, g = droplevels(g)), truncated = FALSE, se = "none"))$loglik)
  dd <- d0; dd$zc <- as.numeric(dd$u)                       # constant within units
  expect_error(cpb_fe(y ~ x + zc, dd, fe = "u"), "do not vary within units")
  expect_error(gec_fe(y ~ x + zc, dd, fe = "u"), "do not vary within units")
})

test_that("missing values are handled by listwise deletion in the fixed-effects and mixture fits", {
  dd <- d0; dd$x[c(3, 30, 100)] <- NA
  fe <- cpb_fe(y ~ x, dd, fe = "u"); expect_equal(fe$n, n - 3L)
  ref <- cpb_fe(y ~ x, dd[!is.na(dd$x), ], fe = "u")
  expect_equal(fe$loglik, ref$loglik, tolerance = 1e-6)
  ge <- gec_fe(y ~ x, dd, fe = "u"); expect_equal(ge$n, n - 3L)
  dz <- d0; dz$z[c(4, 40)] <- NA
  zf <- zi_cpb(y ~ x, dz, zero = ~ z); expect_equal(zf$n, n - 2L)
  zg <- zi_gec(y ~ x, dz, zero = ~ z); expect_equal(zg$n, n - 2L)
})

test_that("predict on the GEC uses the stored factor levels", {
  dg <- transform(d0, g = factor(sample(c("a", "b", "c"), n, TRUE)))
  gg <- gec(y ~ x + g, dg, se = "none")
  p1 <- predict(gg, data.frame(x = 0, g = factor("b", levels = "b")), type = "response")
  p2 <- predict(gg, data.frame(x = 0, g = factor("b", levels = c("a", "b", "c"))), type = "response")
  expect_equal(p1, p2)
  gf <- gec_fe(y ~ x + g, dg, fe = "u")
  expect_equal(predict(gf, data.frame(x = 0, g = factor("b", levels = "b"))),
               predict(gf, data.frame(x = 0, g = factor("b", levels = c("a", "b", "c")))))
})
