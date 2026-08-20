## The over-conditioning guard: (1) a near-saturating mean model SKIPS the
## at-risk screen loudly; (2) a dummy-heavy but fittable design flags the
## verdict as diagnostic (parameter share of positives >= 0.10) and recommends
## the bootstrap threshold; (3) an ordinary design carries no flag.

test_that("saturating design skips the at-risk screen and says why", {
  set.seed(11)
  n <- 40
  d <- data.frame(y = rpois(n, 3) + 1L, g = factor(seq_len(n) %/% 1L))  # one dummy per obs
  s <- suppressWarnings(ud_screen(y ~ g, data = d, run_cpb = FALSE,
                                  run_gp = FALSE, run_comp = FALSE))
  expect_true(s$atrisk_skipped)
  expect_true(s$overconditioned)
  expect_true(is.na(s$pearson_ztp))
  expect_output(print(s), "SKIPPED -- over-conditioning")
})

test_that("dummy-heavy design flags over-conditioning and recommends bootstrap", {
  set.seed(12)
  units <- 30; T_ <- 8; n <- units * T_
  d <- data.frame(y = rpois(n, 4) + 1L, x = rnorm(n),
                  u = factor(rep(seq_len(units), each = T_)))
  s <- suppressWarnings(ud_screen(y ~ x + u, data = d, run_cpb = FALSE,
                                  run_gp = FALSE, run_comp = FALSE))
  ## 31+ fitted params on 240 positives: share > 0.10
  expect_false(s$atrisk_skipped)
  expect_true(is.finite(s$pearson_ztp))
  expect_true(s$overconditioned)
  expect_gte(s$sat_ratio, 0.10)
  expect_output(print(s), "over-conditioning caution")
  expect_output(print(s), "bootstrap")
})

test_that("ordinary design carries no over-conditioning flag", {
  set.seed(13)
  n <- 400
  d <- data.frame(y = rpois(n, 4) + 1L, x = rnorm(n))
  s <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE,
                                  run_gp = FALSE, run_comp = FALSE))
  expect_false(s$atrisk_skipped)
  expect_false(isTRUE(s$overconditioned))
  expect_lt(s$sat_ratio, 0.10)
})
