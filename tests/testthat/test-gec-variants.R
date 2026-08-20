test_that("gec_fe recovers fixed effects and dispersion, and scores", {
  skip_on_cran()
  set.seed(3)
  du <- do.call(rbind, lapply(1:30, function(i) {
    fe <- rnorm(1, 0, 0.5); x <- rnorm(20)
    data.frame(unit = i, x = x, y = rcpb(20, exp(fe + 0.5 * x), 0.5, truncated = FALSE))
  }))
  gf <- gec_fe(y ~ x, du, fe = "unit")
  expect_s3_class(gf, "gec_fe")
  expect_lt(gf$delta, 0.85)                                  # underdispersed panel
  expect_true(is.finite(score(gf)[["logscore"]]))
  expect_length(predict(gf, newdata = data.frame(x = c(-1, 1))), 2L)
})

test_that("hurdle_gec and zi_gec recover their parameters", {
  skip_on_cran()
  set.seed(1); n <- 1000; x <- rnorm(n); z <- rnorm(n)
  dh <- data.frame(y = rhurdle_cpb(n, exp(1.1 + 0.5 * x), 0.5, plogis(-0.2 + 0.8 * z)), x = x, z = z)
  hg <- hurdle_gec(y ~ x, dh, participation = ~ z)
  expect_lt(abs(hg$intensity$coefficients[["x"]] - 0.5), 0.15)
  expect_lt(hg$delta, 0.85)                                  # underdispersed intensity
  expect_true(is.finite(score(hg)[["logscore"]]))

  dz <- data.frame(y = rzicpb(n, exp(1.3 + 0.5 * x), 0.5, plogis(-0.4 + 0.8 * z)), x = x, z = z)
  zg <- zi_gec(y ~ x, dz, zero = ~ z)
  expect_lt(abs(zg$coefficients[["x"]] - 0.5), 0.15)
  expect_lt(zg$delta, 0.85)
  expect_gt(mean(zg$pi_full), 0.2)
})

test_that("zi_test handles the GEC family and enforces the nest", {
  skip_on_cran()
  set.seed(2); n <- 800; x <- rnorm(n)
  dz <- data.frame(y = rzicpb(n, exp(1.3 + 0.4 * x), 0.5, plogis(-0.3)), x = x)
  zg <- zi_gec(y ~ x, dz); g <- gec(y ~ x, dz, se = "none")
  tt <- zi_test(g, zg)
  expect_s3_class(tt, "zi_test")
  expect_identical(tt$family, "GEC")
  ## family mismatch: a CPB is not the nest of a zi_gec
  expect_error(zi_test(cpb(y ~ x, dz, truncated = FALSE, se = "none"), zg), "gec")
  ## two non-inflated models
  expect_error(zi_test(g, gec(y ~ x, dz, se = "none")), "zero-inflated")
})

test_that("compare_models spans CPB and GEC families", {
  skip_on_cran()
  set.seed(4); n <- 600; x <- rnorm(n)
  d <- data.frame(y = rzicpb(n, exp(1.3 + 0.5 * x), 0.5, plogis(-0.3)), x = x)
  cm <- compare_models(gec = gec(y ~ x, d, se = "none"),
                       zi_gec = zi_gec(y ~ x, d),
                       zi_cpb = zi_cpb(y ~ x, d))
  expect_s3_class(cm, "data.frame")
  expect_true(all(c("AIC", "logscore") %in% colnames(cm)))
})
