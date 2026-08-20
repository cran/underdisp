test_that("first_difference decomposes into extensive and intensive margins", {
  set.seed(1); n <- 700; x <- rnorm(n); z <- rnorm(n)
  yh <- rhurdle_cpb(n, exp(1.3 + 0.5 * x), 0.5, plogis(-0.2 + 0.9 * z))
  hf <- hurdle_cpb(y ~ x, data = data.frame(y = yh, x = x, z = z), participation = ~ z)

  # x is an intensity-only covariate: participation change ~ 0, marginal > 0
  fdx <- first_difference(hf, "x", from = -1, to = 1)
  expect_s3_class(fdx, "ud_fd")
  expect_equal(nrow(fdx), 3L)
  expect_named(fdx, c("component", "from", "to", "diff", "lower", "upper", "method"))
  expect_lt(abs(fdx$diff[fdx$component == "participation"]), 1e-8)
  expect_gt(fdx$diff[fdx$component == "marginal"], 0)

  # z is a participation-only covariate: intensity change ~ 0
  fdz <- first_difference(hf, "z", from = -1, to = 1)
  expect_lt(abs(fdz$diff[fdz$component == "intensity"]), 1e-8)
  expect_gt(fdz$diff[fdz$component == "participation"], 0)

  # stage moves the covariate in one equation only
  fs <- first_difference(hf, "z", from = -1, to = 1, stage = "intensity")
  expect_lt(abs(fs$diff[fs$component == "participation"]), 1e-8)
  # unknown arguments error rather than being silently swallowed
  expect_error(first_difference(hf, "z", from = -1, to = 1, stge = "both"),
               "unused argument")

  ii <- irr(hf)
  expect_s3_class(ii, "ud_irr")
  expect_named(ii, c("term", "equation", "ratio", "estimate", "lower", "upper", "method"))
  expect_true(all(c("count", "binary") %in% ii$equation))
  expect_true(all(c("IRR", "OR") %in% ii$ratio))
  ic <- implied_ceiling(hf, newdata = data.frame(x = 0, z = 0))
  expect_true(is.finite(ic$ceiling))
})

test_that("zi_test detects zero-inflation and its absence", {
  skip_on_cran()
  set.seed(2); n <- 600; x <- rnorm(n)
  yz <- rzicpb(n, exp(1.2 + 0.4 * x), 0.5, pi = 0.35)
  tt <- zi_test(y ~ x, data = data.frame(y = yz, x = x))
  expect_s3_class(tt, "zi_test")
  expect_lt(tt$p.value, 0.05)

  yc <- rcpb(n, exp(1.2 + 0.4 * x), 0.5)
  t0 <- zi_test(y ~ x, data = data.frame(y = yc, x = x))
  expect_gt(t0$p.value, 0.05)
})
