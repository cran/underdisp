test_that("cpb supports clustered and ordinary bootstrap", {
  skip_on_cran()
  set.seed(1); n <- 300; x <- rnorm(n); clu <- rep(1:30, each = 10)
  y <- rcpb(n, exp(1 + 0.5 * x), 0.5, truncated = FALSE)
  d <- data.frame(y = y, x = x, clu = clu)
  fc <- cpb(y ~ x, data = d, truncated = FALSE, se = "bootstrap", B = 30, cluster = "clu")
  expect_true(fc$clustered)
  expect_equal(fc$n_clusters, 30L)
  expect_true(is.finite(fc$se.beta[["x"]]))
  fo <- cpb(y ~ x, data = d, truncated = FALSE, se = "bootstrap", B = 30)
  expect_false(fo$clustered)
})

test_that("cpb_fe supports a pairs/cluster bootstrap and is scorable", {
  skip_on_cran()
  set.seed(2)
  d <- do.call(rbind, lapply(1:25, function(u) {
    fe <- rnorm(1, 0, 0.5); x <- rnorm(20)
    data.frame(unit = u, x = x, y = rcpb(20, exp(fe + 0.5 * x), 0.5, truncated = FALSE))
  }))
  f <- cpb_fe(y ~ x, data = d, fe = "unit", se = "bootstrap", B = 15)
  expect_true(is.finite(f$se.beta[["x"]]))
  expect_equal(f$cluster, "unit")
  expect_true(is.finite(score(f)[["logscore"]]))     # cpb_fe now carries its response
})

test_that("zi_cpb supports bootstrap SEs", {
  skip_on_cran()
  set.seed(3); n <- 300; x <- rnorm(n)
  y <- rzicpb(n, exp(1 + 0.4 * x), 0.5, 0.3); d <- data.frame(y = y, x = x)
  f <- zi_cpb(y ~ x, data = d, zero = ~ 1, se = "bootstrap", B = 12)
  expect_true(is.finite(f$se.beta[["x"]]))
  expect_true(is.finite(f$se.zero[[1]]))
})
