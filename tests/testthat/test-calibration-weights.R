## Frequency weights in the calibration tools: a row of weight w counts as w copies of the row, so a weighted fit's
## rootogram, PIT histogram, dispersion profile, and held-out scores are those of the expanded data.
set.seed(31)
n <- 90; x <- seq(-1, 1, length.out = n)
w <- rep(rep_len(c(1, 3, 2), 30), 3)                  # every third of the rows carries the same weight total
d  <- data.frame(y = rcpb(n, exp(1 + 0.8 * x), 0.5, truncated = FALSE), x = x, w = w)
de <- d[rep(seq_len(n), d$w), ]
fw <- cpb(y ~ x, d, truncated = FALSE, se = "none", weights = "w")
fx <- cpb(y ~ x, de, truncated = FALSE, se = "none")
f0 <- cpb(y ~ x, d, truncated = FALSE, se = "none")
f1 <- cpb(y ~ x, d, truncated = FALSE, se = "none", weights = rep(1, n))

test_that("a weighted fit's rootogram and PIT histogram are those of the expanded data", {
  pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off(), add = TRUE)
  rw <- rootogram(fw); rx <- rootogram(fx)
  expect_equal(rw$observed, rx$observed)
  expect_equal(sum(rw$observed), sum(w))
  expect_equal(rw$expected, rx$expected, tolerance = 1e-4)
  expect_equal(pit_hist(fw), pit_hist(fx), tolerance = 1e-4)
  expect_equal(nrow(rootogram(f0, kmax = max(d$y) + 3)), max(d$y) + 4)   # a kmax above the largest count
})

test_that("unit weights reproduce the unweighted tools", {
  pdf(file = tempfile(fileext = ".pdf")); on.exit(grDevices::dev.off(), add = TRUE)
  expect_equal(rootogram(f1), rootogram(f0))
  expect_equal(pit_hist(f1), pit_hist(f0))
  expect_equal(as.data.frame(dispersion_profile(f1, plot = FALSE)), as.data.frame(dispersion_profile(f0, plot = FALSE)))
})

test_that("a weighted dispersion profile holds equal weight per bin and matches the expanded data", {
  expect_false(is.unsorted(fitted(fw)))                # rows in fitted-mean order: the thirds of the weight are the thirds of the rows
  pw <- dispersion_profile(fw, plot = FALSE, bins = 3); px <- dispersion_profile(fx, plot = FALSE, bins = 3)
  expect_equal(pw$n, c(60, 60, 60))
  expect_equal(pw$n, px$n)
  expect_equal(pw$mean_y, px$mean_y)
  expect_equal(pw$mean_fitted, px$mean_fitted, tolerance = 1e-4)
  expect_equal(pw$ratio_empirical, px$ratio_empirical, tolerance = 1e-4)
  expect_equal(pw$ratio_cpb, px$ratio_cpb, tolerance = 1e-4)
})

test_that("held-out scores weight the rows of newdata and the folds of cv_score", {
  expect_equal(score(fw, newdata = d, weights = "w"), score(fw), tolerance = 1e-8)
  expect_error(score(fw, weights = "w"), "newdata")
  skip_on_cran()
  folds <- rep_len(1:3, n)
  cw <- cv_score(function(tr) cpb(y ~ x, tr, truncated = FALSE, se = "none", weights = "w"), d, folds = folds, weights = "w")
  cx <- cv_score(function(tr) cpb(y ~ x, tr, truncated = FALSE, se = "none"), de, folds = rep(folds, d$w))
  expect_equal(cw$logscore, cx$logscore, tolerance = 1e-4)
  expect_equal(cw$rps, cx$rps, tolerance = 1e-4)
})
