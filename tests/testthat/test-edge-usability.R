## Usability at the edges: unused factor levels are dropped as lm() and glm() drop them, an all-zero response
## warns that the dispersion parameter is not identified, and newdata rows with a missing covariate predict NA.
set.seed(11)
n <- 300; x <- rnorm(n); g <- factor(sample(letters[1:3], n, TRUE))
d <- data.frame(y = rcpb(n, exp(1 + 0.4 * x), 0.5), x = x, g = g)
du <- d; du$g <- factor(du$g, levels = c(levels(d$g), "zz"))

test_that("unused factor levels are dropped, as lm() and glm() drop them", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  a <- cpb(y ~ x + g, du, truncated = FALSE, se = "none")
  b <- cpb(y ~ x + g, d, truncated = FALSE, se = "none")
  expect_equal(a$loglik, b$loglik)
  expect_equal(coef(a), coef(b))
  expect_equal(predict(a, newdata = d[1:5, ]), predict(b, newdata = d[1:5, ]))
  expect_equal(count_reg(y ~ x + g, du, se = "none")$loglik, count_reg(y ~ x + g, d, se = "none")$loglik)
})

test_that("an all-zero response fits with a warning that the dispersion is not identified", {
  z0 <- data.frame(y = rep(0L, 40), x = rnorm(40))
  expect_warning(cpb(y ~ x, z0, truncated = FALSE, se = "none"), "every count is zero")
  expect_warning(count_reg(y ~ x, z0, se = "none"), "every count is zero")
})

test_that("newdata rows with a missing covariate predict NA", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  nd <- data.frame(x = c(0, NA, 1), g = factor(c("a", "b", "c"), levels = levels(d$g)))
  fits <- list(cpb(y ~ x + g, d, truncated = FALSE, se = "none"),
               count_reg(y ~ x + g, d, family = "compois", se = "none"),
               count_reg(y ~ x + g, d, family = "poisson", se = "none"))
  for (f in fits) {
    p <- predict(f, newdata = nd)
    expect_true(is.na(p[2]))
    expect_true(all(is.finite(p[c(1, 3)])))
  }
  f <- fits[[1]]
  expect_true(is.na(predict(f, newdata = data.frame(x = NA, g = "a"))))
  expect_error(predict(f, newdata = data.frame(x = TRUE, g = "a")), "type it had in the fitting data")
})

test_that("fixed-effects and zero-inflated predictions also return NA for a missing covariate", {
  skip_on_cran()
  dd <- transform(d, u = factor(rep_len(1:15, n)), z = rnorm(n))
  ff <- cpb_fe(y ~ x, dd, fe = "u", se = "none")
  p <- predict(ff, newdata = data.frame(x = c(0, NA)))
  expect_true(is.finite(p[1]) && is.na(p[2]))
  zc <- zi_count(y ~ x, dd, family = "poisson", zero = ~ z)
  pz <- predict(zc, newdata = data.frame(x = c(0, NA), z = c(0, 0)))
  expect_true(is.finite(pz[1]) && is.na(pz[2]))
})
