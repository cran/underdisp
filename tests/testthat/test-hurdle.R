test_that("hurdle_cpb recovers participation and intensity parameters", {
  set.seed(1)
  n <- 1200; x <- rnorm(n); z <- rnorm(n)
  y <- rhurdle_cpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5,
                   p = plogis(-0.2 + 0.8 * z))
  fit <- hurdle_cpb(y ~ x, data = data.frame(y = y, x = x, z = z),
                    participation = ~ z, se = "none")
  expect_s3_class(fit, "hurdle_cpb")

  # participation logit recovers the sign and rough magnitude of z
  pc <- coef(fit$participation)
  expect_true(pc[["z"]] > 0.4 && pc[["z"]] < 1.2)

  # intensity recovers the slope and an underdispersed alpha
  expect_true(fit$intensity$coefficients[["x"]] > 0.3 &&
              fit$intensity$coefficients[["x"]] < 0.7)
  expect_true(fit$intensity$alpha > 0 && fit$intensity$alpha < 1)

  # predictions decompose: marginal = participation * intensity
  nd <- data.frame(x = 0, z = 0)
  p  <- predict(fit, nd, type = "participation")
  ey <- predict(fit, nd, type = "intensity")
  m  <- predict(fit, nd, type = "response")
  expect_equal(unname(m), unname(p * ey), tolerance = 1e-8)
  expect_true(p > 0 && p < 1 && ey >= 1)
})

test_that("rcpb and rhurdle_cpb behave", {
  set.seed(2)
  expect_true(all(rcpb(200, 3, 0.5, truncated = TRUE) >= 1))
  expect_error(rcpb(10, 3, 1.5))            # alpha out of range
  y <- rhurdle_cpb(500, lambda = 2, alpha = 0.4, p = 0.5)
  expect_true(mean(y == 0) > 0.2)           # a genuine hurdle at zero
})

test_that("the zero-truncated CPB at a ceiling below 1 is the point mass at 1 everywhere", {
  ## lambda = 0.3, alpha = 0.5: ceiling 0.6, so the untruncated CPB is the point mass at 0
  expect_equal(underdisp:::.cpb_pmf1(0.3, 0.5, 3, truncated = TRUE), c(0, 1, 0, 0))
  expect_equal(dcpb(0:2, 0.3, 0.5, truncated = TRUE), c(0, 1, 0))
  expect_equal(underdisp:::.cpb_ztmean(0.3, 0.5), 1)
  expect_equal(underdisp:::.cpb_mean1(0.3, 0.5, truncated = TRUE), 1)
  expect_equal(unname(underdisp:::.cpb_moments(0.3, 0.5, truncated = TRUE)$var), 0)
  ## the simulator follows the same rule and says so: the stated parameters give such a count probability zero
  set.seed(1); expect_warning(yy <- rcpb(20, 0.3, 0.5, truncated = TRUE), "ceiling")
  expect_true(all(yy == 1))
  expect_silent(rcpb(20, 3, 0.5, truncated = TRUE))
  ## at a ceiling of exactly 1 the same distribution holds
  expect_equal(underdisp:::.cpb_pmf1(0.5, 0.5, 3, truncated = TRUE), c(0, 1, 0, 0))
})
