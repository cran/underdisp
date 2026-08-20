test_that("gec recovers the dispersion direction across the range", {
  skip_on_cran()
  set.seed(1); n <- 1200; x <- rnorm(n)
  gp <- gec(y ~ x, data.frame(y = rpois(n, exp(1 + 0.5 * x)), x = x), se = "none")
  gn <- gec(y ~ x, data.frame(y = rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)), x = x), se = "none")
  gu <- gec(y ~ x, data.frame(y = rcpb(n, exp(1.3 + 0.5 * x), 0.5, truncated = FALSE), x = x), se = "none")
  expect_gt(gp$delta, 0.85); expect_lt(gp$delta, 1.20)   # Poisson -> ~ equidispersed
  expect_gt(gn$delta, 1.30)                              # NB -> overdispersed
  expect_lt(gu$delta, 0.70)                              # CPB -> underdispersed
  expect_true(all(abs(c(gp$coefficients[["x"]], gn$coefficients[["x"]], gu$coefficients[["x"]]) - 0.5) < 0.12))
})

test_that("gec integrates with compare_dispersion, score, broom, compare_models", {
  skip_on_cran()
  set.seed(2); n <- 600; x <- rnorm(n)
  d <- data.frame(y = rcpb(n, exp(1.3 + 0.5 * x), 0.5, truncated = FALSE), x = x)
  cmp <- compare_dispersion(y ~ x, d)
  expect_true("GEC" %in% rownames(cmp$table))
  expect_true(is.finite(cmp$gec_delta))
  g <- gec(y ~ x, d, se = "none")
  expect_true(is.finite(score(g)[["logscore"]]))
  cm <- compare_models(gec = g, cpb = cpb(y ~ x, d, truncated = FALSE, se = "none"))
  expect_s3_class(cm, "data.frame")
})

test_that("gec validates newdata and predicts", {
  skip_on_cran()
  set.seed(3); n <- 400; x <- rnorm(n)
  g <- gec(y ~ x, data.frame(y = rpois(n, exp(1 + 0.5 * x)), x = x), se = "none")
  expect_error(predict(g, newdata = data.frame(w = 1)), "missing required")
  expect_length(predict(g, newdata = data.frame(x = c(-1, 1))), 2L)
  expect_error(gec(I(x) ~ x, data.frame(x = x), se = "none"), "non-negative integer")  # non-count response
})
