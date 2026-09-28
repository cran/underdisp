## The support guard's warning comes once, from the fit the user asked for, and not from the internal fits of a
## test, a start value, or a bootstrap or jackknife replicate.
collect <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, warnings = msgs)
}

test_that("the guard filter muffles only the guard warning", {
  expect_warning(underdisp:::.ud_quiet_guard(warning("something else")), "something else")
  expect_silent(underdisp:::.ud_quiet_guard(warning("the fitted ceiling reaches max.support (10): alpha is bounded")))
})

test_that("zi_test() from a formula does not pass on its nest's guard warning", {
  skip_on_cran()
  set.seed(1); x <- rnorm(500)
  y <- rzicpb(500, lambda = exp(1.2 + 0.4 * x), alpha = 0.5, pi = 0.3)
  r <- collect(zi_test(y ~ x, data = data.frame(y = y, x = x)))
  expect_s3_class(r$value, "zi_test")
  expect_false(any(grepl("max.support", r$warnings, fixed = TRUE)))
})

test_that("a cpb_fe bootstrap reports the guard once, not once per replicate", {
  skip_on_cran()
  set.seed(7); x <- rnorm(120)                            # overdispersed counts: the ceiling runs to the guard
  d <- data.frame(y = rnbinom(120, mu = exp(1 + 0.3 * x), size = 1), x = x, u = rep(1:12, each = 10))
  r <- collect(cpb_fe(y ~ x, d, fe = "u", se = "bootstrap", B = 4, max.support = 60))
  expect_true(r$value$support_binding)
  expect_identical(sum(grepl("reaches max.support", r$warnings, fixed = TRUE)), 1L)
})
