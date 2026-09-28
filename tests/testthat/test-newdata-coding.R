## predict()/score() on newdata reproduce the fit's factor coding, including a
## reference level set through the factor's contrasts attribute (which is lost
## when model.frame() re-levels a factor through xlev; AER's NMES1988 carries
## such a factor).

set.seed(19); n <- 600; x <- rnorm(n)
g <- factor(sample(c("poor", "average", "excellent"), n, TRUE), levels = c("poor", "average", "excellent"))
contrasts(g) <- stats::contr.treatment(3, base = 2)          # reference level "average"
z <- rnorm(n)
d <- data.frame(y = rcpb(n, exp(0.8 + 0.4 * x + 0.3 * (g == "excellent")), 0.5), x = x, g = g, z = z)
d$y0 <- ifelse(runif(n) < plogis(-1 + 0.5 * z), 0L, d$y)
sub <- d[1:80, ]

test_that("single-equation fits keep the coding on newdata", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  for (f in list(cpb(y ~ x + g, d, truncated = FALSE, se = "none"),
                 count_reg(y ~ x + g, d, family = "poisson"),
                 gec(y ~ x + g, d, se = "none"))) {
    expect_true(any(grepl("^gpoor$|^g1$", names(coef(f)))))     # the coding with "average" as reference
    expect_equal(unname(predict(f, newdata = sub)), unname(fitted(f)[1:80]), tolerance = 1e-10)
    expect_true(all(is.finite(score(f, newdata = sub))))
  }
})

test_that("two-part fits keep the coding on newdata in both equations", {
  skip_on_cran()                                  # heavy: runs in CI and locally, not on CRAN's clock
  z1 <- zi_cpb(y0 ~ x, d, zero = ~ g, se = "none")
  expect_equal(unname(predict(z1, newdata = d)), unname(fitted(z1)), tolerance = 1e-10)
  z2 <- zi_count(y0 ~ x + g, d, family = "poisson", zero = ~ g)
  expect_equal(unname(predict(z2, newdata = d)), unname(fitted(z2)), tolerance = 1e-10)
  h <- hurdle_cpb(y0 ~ x + g, d, participation = ~ g, se = "none")
  expect_equal(unname(predict(h, newdata = d)), unname(fitted(h)), tolerance = 1e-10)
  expect_true(all(is.finite(score(z1, newdata = sub)))); expect_true(all(is.finite(score(h, newdata = sub))))
})

test_that("a newdata factor with the wrong levels is refused, not multiplied by position", {
  f <- count_reg(y ~ x + g, d, family = "poisson")
  bad <- sub; bad$g <- factor(as.character(bad$g), levels = c("excellent", "poor", "average"))
  ## the same levels in a different order still code correctly through xlev
  expect_equal(unname(predict(f, newdata = bad)), unname(fitted(f)[1:80]), tolerance = 1e-10)
  bad$g <- factor(sample(c("a", "b"), 80, TRUE))
  expect_error(predict(f, newdata = bad))
})
