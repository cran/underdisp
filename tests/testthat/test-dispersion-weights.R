test_that("the parametric bootstrap treats a weight of w as w independent observations", {
  skip_on_cran()
  set.seed(71); n <- 50; x <- rnorm(n)
  d <- data.frame(y = rpois(n, exp(0.8 + 0.3 * x)), x = x); d$w <- 2
  dd <- d[rep(seq_len(n), each = 2), ]                           # the same rows, written out, in the expansion's order
  fw <- count_reg(y ~ x, d, family = "genpois", weights = "w")
  fe <- count_reg(y ~ x, dd, family = "genpois")
  set.seed(9); tw <- dispersion_test(fw, B = 4)       # four replicates are enough to compare draw for draw
  set.seed(9); te <- dispersion_test(fe, B = 4)
  expect_equal(unname(tw$statistic), unname(te$statistic), tolerance = 1e-6)
  ## same seed, same expanded rows: the replicate statistics agree draw for draw (one response per weighted row,
  ## with the weight passed to the refit, doubled them on average)
  expect_equal(tw$boot, te$boot, tolerance = 1e-5)
  expect_equal(tw$p.value, te$p.value)
})

test_that("weights that are not whole numbers skip the bootstrap with a note", {
  set.seed(72); n <- 80; x <- rnorm(n)
  d <- data.frame(y = rpois(n, exp(0.8 + 0.3 * x)), x = x); d$w <- runif(n, 0.5, 1.5)
  f <- count_reg(y ~ x, d, family = "genpois", weights = "w")
  t1 <- dispersion_test(f, B = 29)
  expect_identical(t1$calibration, "asymptotic")
  expect_true(any(grepl("not whole numbers", t1$note)))
  expect_true(is.finite(t1$p.value))
})
