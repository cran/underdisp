## The COM-Poisson is parameterized on its rate. On a nearly degenerate series its likelihood rises along a
## ridge on which nu grows and the rate, about the mode to the power nu, grows without bound, so the fit has
## to reach the likelihood's supremum there instead of stopping where the rate crosses a fixed bound.
saturated_ll <- function(y) { tb <- table(y); sum(tb * log(tb / length(y))) }

test_that("the excursion bound applies on the COM-Poisson's mean scale", {
  cmp <- underdisp:::.count_fam("compois"); poi <- underdisp:::.count_fam("poisson")
  expect_false(underdisp:::.count_excursion(118, 74, cmp))       # rate e^118, mean about e^1.6
  expect_true(underdisp:::.count_excursion(20, 1, cmp))          # mean e^20
  expect_true(underdisp:::.count_excursion(20, NA_real_, poi))
  expect_false(underdisp:::.count_excursion(log(5), NA_real_, poi))
})

test_that("the zero-truncated COM-Poisson reaches the saturated likelihood on two-valued series", {
  skip_on_cran()
  for (y in list(c(1L, rep(2L, 63)), c(rep(4L, 47), rep(5L, 16)))) {
    m <- count_reg(y ~ 1, data = data.frame(y = y), family = "compois", truncated = TRUE, se = "none")
    expect_true(m$converged)
    expect_equal(m$loglik, saturated_ll(y), tolerance = 1e-5)
  }
  ## an interior optimum is unaffected
  y <- rep(1:3, c(41, 109, 45))
  m <- count_reg(y ~ 1, data = data.frame(y = y), family = "compois", truncated = TRUE, se = "none")
  expect_equal(m$loglik, -196.8499, tolerance = 1e-6)
})
