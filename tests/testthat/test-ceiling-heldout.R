## The CPB's bound as a quantity (implied_ceiling()'s max_count) and as a held-out test (ceiling_exceedance());
## first_difference() bootstraps that refit with the fit's own settings; the COM-Poisson and hyper-Poisson shape
## parameters held at the family's limits.

test_that("implied_ceiling() reports the largest feasible count on each support", {
  skip_on_cran()
  set.seed(41); x <- rnorm(300); d <- data.frame(y = rcpb(300, exp(1.2 + 0.4 * x), 0.5), x = x)
  for (s in c("ceiling", "floor")) {
    f <- cpb(y ~ x, d, truncated = FALSE, se = "none", support = s)
    ic <- implied_ceiling(f)
    expect_identical(ic$max_count, if (s == "ceiling") ceiling(ic$ceiling - 1e-9) else floor(ic$ceiling + 1e-9))
    expect_length(f$Y, 300L)
    expect_true(all(f$Y <= ic$max_count))                       # in sample, every feasible fit respects its bound
  }
  z <- zi_cpb(y ~ x, d, zero = ~ 1)
  expect_equal(implied_ceiling(z)$ceiling, z$lambda_full / (1 - z$alpha))   # the rate, not the mean, on King's support
  expect_equal(implied_ceiling(z, newdata = d[1:5, ])$lambda, predict(z, newdata = d[1:5, ], type = "rate"))
})

test_that("ceiling_exceedance() counts held-out counts above the fitted bound, at each unit's own effect", {
  set.seed(42); x <- rnorm(400); d <- data.frame(y = rcpb(400, exp(1.2 + 0.4 * x), 0.5), x = x)
  f <- cpb(y ~ x, d[1:300, ], truncated = FALSE, se = "none")
  ho <- d[301:400, ]
  top <- implied_ceiling(f, newdata = ho)$max_count
  ho$y[1:3] <- top[1:3] + 1L                                    # three counts the fitted CPB cannot produce
  ce <- ceiling_exceedance(f, ho)
  expect_equal(ce$exceed, 3 + sum(ho$y[-(1:3)] > top[-(1:3)]))
  expect_equal(ce$share, ce$exceed / 100)
  expect_equal(ceiling_exceedance(f, ho, weights = rep(2, 100))$share, ce$share)
  ## fixed effects: a unit seen in training enters at its own effect
  set.seed(43); u <- rep(1:20, each = 15); xx <- rnorm(300); ue <- rnorm(20, 0, 0.8)
  p <- data.frame(y = rcpb(300, exp(1 + ue[u] + 0.3 * xx), 0.5), x = xx, u = u, t = rep(1:15, 20))
  fe <- cpb_fe(y ~ x, p[p$t <= 10, ], fe = "u")
  hp <- p[p$t > 10, ]
  own <- fe$fe[as.character(hp$u)]
  top <- ceiling(exp(own + fe$coefficients[["x"]] * hp$x) / (1 - fe$alpha) - 1e-9)
  expect_equal(ceiling_exceedance(fe, hp)$exceed, sum(hp$y > top))
  expect_error(ceiling_exceedance(count_reg(y ~ x, d, family = "poisson"), ho), "CPB family")
})

test_that("compare_dispersion() and ud_screen() no longer report an in-sample exceedance share", {
  set.seed(44); x <- rnorm(200); d <- data.frame(y = rcpb(200, exp(1.4 + 0.3 * x), 0.5), x = x)
  cd <- compare_dispersion(y ~ x, d)
  expect_null(cd$ceiling_exceedance)
  sc <- ud_screen(y ~ x, data = d)
  expect_null(sc$ceil_exceed)
  expect_output(print(sc), "CPB \\(full data\\): alpha")
})

test_that("first_difference() bootstraps refit with the fit's fixed effects, link, and offset", {
  skip_on_cran()
  set.seed(45); n <- 400; u <- rep(1:20, each = 20); x <- rnorm(n); z <- rnorm(n); ex <- runif(n, 0.5, 2)
  y <- ifelse(rbinom(n, 1, plogis(0.2 + 0.7 * z)) == 1,
              rcpb(n, ex * exp(0.8 + rnorm(20, 0, 0.5)[u] + 0.4 * x), 0.5, truncated = TRUE), 0L)
  d <- data.frame(y = y, x = x, z = z, u = u, lex = log(ex))
  h <- hurdle_cpb(y ~ x, d, participation = ~ z, fe = "u", link = "probit", offset = "lex")
  ra <- underdisp:::.fd_refit_args(h, c("fe", "part_fe", "link"))
  expect_identical(ra$fe, "u"); expect_identical(ra$link, "probit"); expect_identical(ra$offset, "lex")
  set.seed(1); fd <- first_difference(h, "x", from = -1, to = 1, B = 4, data = d)
  expect_s3_class(fd, "ud_fd"); expect_true(all(is.finite(fd$lower)))
  z <- suppressWarnings(zi_cpb(y ~ x, d, zero = ~ z, fe = "u", zero_fe = "u", method = "em"))   # settings, not the fit
  expect_identical(z$refit_args$zero_fe, "u"); expect_identical(z$refit_args$method, "em")
  ## a vector offset cannot follow resampled rows
  h2 <- hurdle_cpb(y ~ x, d, participation = ~ z, offset = d$lex)
  expect_error(first_difference(h2, "x", from = -1, to = 1, B = 2, data = d), "column name")
  ## an object saved before 0.2.0 (no refit_args) still refits, with the fields it carries
  h3 <- h; h3$refit_args <- NULL
  expect_s3_class(first_difference(h3, "x", from = -1, to = 1, B = 2, data = d), "ud_fd")
})

test_that("the COM-Poisson and hyper-Poisson shapes stop at the family's limits, with a note", {
  skip_on_cran()
  set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1, rgammacount(n, exp(1 + 0.4 * x), 2) + 1L, 0L)
  d <- data.frame(y = y, x = x)
  for (fam in c("mpcmp", "hyperpois")) {
    t <- system.time(m <- count_reg(y ~ x, d, family = fam))[["elapsed"]]
    expect_true(m$theta >= 1e-8 && m$theta <= 1e8)
    expect_output(print(m), "\\[at the geometric limit\\]")
    expect_lt(t, 300)
  }
  expect_identical(underdisp:::.shape_inv_bounded(c(-50, 0, 50)), c(1e-8, 1, 1e8))
})
