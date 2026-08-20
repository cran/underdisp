test_that("broom methods cover the whole family", {
  skip_on_cran(); skip_if_not_installed("broom")
  set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
  y <- rhurdle_cpb(n, exp(1 + 0.5 * x), 0.5, plogis(-0.2 + 0.8 * z))
  d <- data.frame(y = y, x = x, z = z)
  h  <- hurdle_cpb(y ~ x, data = d, participation = ~ z, se = "none")
  zf <- zi_cpb(y ~ x, data = d, zero = ~ z)
  th <- broom::tidy(h)
  expect_true(all(c("component", "term", "estimate") %in% names(th)))
  expect_true(all(c("participation", "intensity") %in% th$component))
  expect_true(all(c("count", "zero") %in% broom::tidy(zf)$component))
  expect_s3_class(broom::glance(h), "data.frame")
})

test_that("texreg extract methods build a texreg object", {
  skip_on_cran(); skip_if_not_installed("texreg")
  set.seed(2); n <- 300; x <- rnorm(n)
  cf <- cpb(y ~ x, data = data.frame(y = rcpb(n, exp(1 + 0.5 * x), 0.5, truncated = FALSE), x = x),
            truncated = FALSE, se = "none")
  tr <- texreg::extract(cf)
  expect_s4_class(tr, "texreg")
  expect_length(tr@coef, 2L)
})
