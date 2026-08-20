test_that("hurdle_cpb absorbs intensity fixed effects", {
  skip_on_cran()
  set.seed(1)
  d <- do.call(rbind, lapply(1:60, function(u) {
    fe_u <- rnorm(1, 0, 0.6); x <- rnorm(25); z <- rnorm(25)
    data.frame(unit = u, x = x, z = z,
               y = rhurdle_cpb(25, exp(fe_u + 0.5 * x), 0.5, plogis(-0.3 + 0.9 * z)))
  }))
  hf <- hurdle_cpb(y ~ x, data = d, participation = ~ z, fe = "unit")
  expect_s3_class(hf, "hurdle_cpb")
  expect_s3_class(hf$intensity, "cpb_fe")
  expect_true(hf$intensity$alpha > 0 && hf$intensity$alpha < 1)     # underdispersed
  expect_true(coef(hf$participation)[["z"]] > 0.4)                  # participation signal recovered

  fd <- first_difference(hf, "x", from = -1, to = 1)               # decomposition still works under FE
  expect_equal(nrow(fd), 3L)
  expect_true(all(is.finite(fd$diff)))
})

test_that("zi_cpb absorbs intensity fixed effects via classification EM", {
  skip_on_cran()
  set.seed(2)
  d <- do.call(rbind, lapply(1:60, function(u) {
    fe_u <- rnorm(1, 0, 0.6); x <- rnorm(25); z <- rnorm(25)
    data.frame(unit = u, x = x, z = z,
               y = rzicpb(25, exp(fe_u + 0.5 * x), 0.5, plogis(-0.4 + 0.7 * z)))
  }))
  zf <- zi_cpb(y ~ x, data = d, zero = ~ z, fe = "unit")
  expect_s3_class(zf, "zi_cpb")
  expect_false(is.null(zf$fe_hat))
  expect_true(zf$alpha > 0 && zf$alpha < 1)
  expect_true(zf$zero_coef[["z"]] > 0.3)
  expect_length(zf$lambda_full, nrow(d))
})
