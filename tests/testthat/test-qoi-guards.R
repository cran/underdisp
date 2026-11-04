## Quantities of interest that would otherwise be silently wrong: the hurdle decomposition under non-logit links,
## first differences of a variable that enters other columns, ratio labels under non-logit links, the alpha
## interval of a fixed-effects fit, and the two-part decomposition's documentation of its rows.

set.seed(3); n_g <- 400; x_g <- rnorm(n_g); z_g <- rnorm(n_g)
d_g <- data.frame(y = rhurdle_cpb(n_g, lambda = exp(1 + 0.3 * x_g), alpha = 0.5, p = pnorm(0.2 + 0.8 * z_g)),
                  x = x_g, z = z_g)

test_that("hurdle_cpb's first difference uses the participation link it was fitted with", {
  for (lk in c("probit", "cloglog")) {
    h <- hurdle_cpb(y ~ x, d_g, participation = ~ z, link = lk)
    fd <- first_difference(h, "z", 0, 1)
    pp <- predict(h, data.frame(x = mean(d_g$x), z = c(0, 1)), type = "participation")
    expect_equal(c(fd$from[1], fd$to[1]), pp, tolerance = 1e-8)
  }
})

test_that("first_difference() refuses a variable that also enters a transformed or interaction column", {
  skip_on_cran()
  cr <- count_reg(y ~ x + I(x^2), d_g, family = "poisson")
  expect_error(first_difference(cr, "x", 0, 2), "also enters the term")
  cl <- count_reg(y ~ x + log(x + 5), d_g, family = "poisson")
  expect_error(first_difference(cl, "x", 0, 2), "log")
  ci <- count_reg(y ~ x * z, d_g, family = "poisson")
  expect_error(first_difference(ci, "x", 0, 1), "x:z")
  ok <- count_reg(y ~ x + z, d_g, family = "poisson")
  expect_s3_class(first_difference(ok, "x", 0, 1), "ud_fd")
  zg <- zi_gec(y ~ x * z, d_g, zero = ~ z)
  expect_error(first_difference(zg, "x", 0, 1), "x:z")
})

test_that("irr() reports odds ratios only under the logit link", {
  hp <- hurdle_count(y ~ x, d_g, family = "poisson", participation = ~ z, link = "probit")
  expect_message(r <- irr(hp), "not an odds ratio")
  expect_false("OR" %in% r$ratio)
  hl <- hurdle_count(y ~ x, d_g, family = "poisson", participation = ~ z)
  expect_true("OR" %in% irr(hl)$ratio)
})

test_that("alpha_confint() refuses a fixed-effects fit", {
  set.seed(4)
  pan <- do.call(rbind, lapply(1:20, function(i) { xx <- rnorm(8); data.frame(u = i, x = xx, y = rcpb(8, exp(1 + .3 * xx), .5)) }))
  expect_error(alpha_confint(cpb_fe(y ~ x, pan, fe = "u")), "takes a \"cpb\" fit")
})

test_that("the GEC's quantile and random functions refuse where the support guard binds, and infinite rates give NA", {
  expect_warning(r <- rgec(5, lambda = 1000, delta = 1), "binds")
  expect_true(all(is.na(r)))
  expect_warning(q <- qgec(0.5, 1000, 1), "binds")
  expect_true(is.na(q))
  expect_identical(qgec(c(0.1, 0.5, 0.9), 3, 0.7), c(1, 3, 5))
  m <- underdisp:::gec_mean_cpp(c(1, Inf, NA, 0), 0.5, 500L)
  expect_true(is.finite(m[1]) && is.na(m[2]) && is.na(m[3]) && m[4] == 0)
})

test_that("King's CPB pmf in R is refused, not renormalized, past its cap", {
  expect_warning(d <- dcpb(2e5, 2e5, 0.5), "1e5-term cap")
  expect_true(is.na(d))
  set.seed(9); xx <- rnorm(300); d3 <- data.frame(y = rcpb(300, exp(log(400) + .1 * xx), .5), x = xx)
  f <- cpb(y ~ x, d3, truncated = FALSE, max.support = 300, se = "none")
  expect_gt(max(simulate(f, nsim = 1, seed = 1)[[1]]), 300)            # max.support no longer cuts King's support
})
