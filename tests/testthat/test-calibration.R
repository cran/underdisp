test_that("score, rootogram, and pit_hist work on a CPB fit", {
  set.seed(1)
  y <- rcpb(600, lambda = 3, alpha = 0.5)
  fit <- cpb(y ~ 1, data = data.frame(y = y), truncated = FALSE, se = "none")

  s <- score(fit)
  expect_named(s, c("logscore", "rps"))
  expect_true(all(s > 0 & is.finite(s)))

  # direct all plotting to a temporary device
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  rg <- rootogram(fit)
  expect_true(all(c("count", "observed", "expected") %in% names(rg)))
  expect_equal(sum(rg$observed), length(y))
  # expected total is close to observed total (mass is nearly all captured)
  expect_true(abs(sum(rg$expected) - length(y)) < 0.05 * length(y))

  ph <- pit_hist(fit, bins = 10)
  expect_length(ph, 10)
  expect_true(all(ph >= 0))
})

test_that("score prefers the correct model to a misspecified Poisson", {
  set.seed(3)
  y <- rcpb(800, lambda = 3, alpha = 0.5)
  d <- data.frame(y = y)
  cpb_fit <- cpb(y ~ 1, data = d, truncated = FALSE, se = "none")
  # a crude Poisson log score on the same data
  lam <- mean(y); kmax <- max(y)
  P <- matrix(dpois(0:kmax, lam), length(y), kmax + 1, byrow = TRUE)
  ls_pois <- -mean(log(pmax(P[cbind(seq_along(y), pmin(y, kmax) + 1)], 1e-12)))
  expect_lt(score(cpb_fit)[["logscore"]], ls_pois)
})
