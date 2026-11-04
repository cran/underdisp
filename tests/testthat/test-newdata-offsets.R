## A fit made with an offset carries it to new rows: predict(), score() and ceiling_exceedance() read the fit's
## offset column from newdata (or take `offset =`), refuse rather than assume exposure 1, and first differences
## hold the offset at its mean in every class.

set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n); lexpo <- log(runif(n, 1, 4))
d_off <- data.frame(y = rcpb(n, exp(0.8 + 0.3 * x + lexpo), 0.5), x = x, z = z, lexpo = lexpo)

test_that("new rows take the fit's offset column, and a fit with an offset refuses rows without one", {
  f <- cpb(y ~ x, d_off[1:300, ], truncated = FALSE, offset = "lexpo", se = "none")
  expect_equal(predict(f, newdata = d_off[1:300, ], type = "rate"), f$rate)
  expect_identical(ceiling_exceedance(f, d_off[301:400, ])$exceed, 0)          # the model's own draws respect the bound
  expect_error(predict(f, newdata = d_off[1:5, c("y", "x")]), "fitted with an offset")
  expect_equal(predict(f, newdata = d_off[1:5, c("y", "x")], offset = d_off$lexpo[1:5], type = "rate"), f$rate[1:5])
  fv <- cpb(y ~ x, d_off, truncated = FALSE, offset = d_off$lexpo, se = "none")  # a vector offset has no column
  expect_error(predict(fv, newdata = d_off[1:5, ]), "fitted with an offset")
  sc <- score(f, newdata = d_off[301:400, ])
  expect_true(is.finite(sc[["logscore"]]))
  g <- gec(y ~ x, d_off, offset = "lexpo", se = "none")
  expect_equal(predict(g, newdata = d_off, type = "link"), g$linear.predictors)
})

test_that("zi_cpb carries its offset to new rows and to the bound check", {
  zo <- suppressWarnings(zi_cpb(y ~ x, d_off, zero = ~ 1, offset = "lexpo"))
  expect_equal(predict(zo, newdata = d_off, type = "rate"), zo$lambda_full)
  expect_identical(ceiling_exceedance(zo, d_off)$exceed, 0)
  expect_equal(zo$int_fe_ref, mean(d_off$lexpo))
})

test_that("two-part first differences hold the offset at its mean", {
  d <- d_off; d$y0 <- ifelse(runif(n) < 0.3, 0L, d$y)
  zc <- zi_count(y0 ~ x, d, family = "poisson", zero = ~ z, offset = "lexpo")
  fd <- first_difference(zc, "x", 0, 1)
  expect_equal(fd$from[fd$component == "count"], unname(exp(mean(d$lexpo) + zc$coefficients[1])))
  hc <- hurdle_count(y0 ~ x, d, family = "poisson", participation = ~ z, offset = "lexpo")
  expect_equal(hc$int_fe_ref, mean(d$lexpo))
  hg <- hurdle_gec(y0 ~ x, d, participation = ~ z, offset = "lexpo")
  fh <- first_difference(hg, "x", 0, 0)
  expect_equal(fh$from[fh$component == "intensity"],
               predict(hg, newdata = data.frame(x = 0, z = mean(d$z), lexpo = mean(d$lexpo)), type = "intensity"),
               tolerance = 1e-6)
})

test_that("score() and cv_score() take the held-out rows' offset", {
  skip_on_cran()
  fc <- cpb(y ~ x, d_off[1:300, ], truncated = FALSE, offset = "lexpo", se = "none")
  fv <- cpb(y ~ x, d_off[1:300, ], truncated = FALSE, offset = d_off$lexpo[1:300], se = "none")
  expect_error(score(fv, newdata = d_off[301:400, ]), "fitted with an offset")
  expect_equal(score(fv, newdata = d_off[301:400, ], offset = d_off$lexpo[301:400]),
               score(fc, newdata = d_off[301:400, ]))
  expect_error(score(fc, offset = "lexpo"), "rows of 'newdata'")
  folds <- rep(1:4, length.out = n)
  cv_c <- cv_score(function(tr) cpb(y ~ x, tr, truncated = FALSE, offset = "lexpo", se = "none"), d_off,
                   kmax = 40, folds = folds)
  vfit <- function(tr) cpb(y ~ x, tr, truncated = FALSE, offset = tr$lexpo, se = "none")
  cv_v <- cv_score(vfit, d_off, kmax = 40, folds = folds, offset = d_off$lexpo)
  expect_equal(cv_v$n_scored, n)
  expect_equal(cv_v$logscore, cv_c$logscore)
  msg <- character(0)
  cv_0 <- withCallingHandlers(cv_score(vfit, d_off, kmax = 40, folds = folds),
                              warning = function(w) { msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_identical(cv_0$n_scored, 0L)
  expect_true(length(msg) == 4L && all(grepl("fitted with an offset", msg)))
})

test_that("a unit a fixed-effects hurdle has not seen takes its own offset once", {
  set.seed(12); Tt <- 6
  mk <- function(units) {
    m <- length(units) * Tt
    d <- data.frame(unit = rep(units, each = Tt), x = rnorm(m), lexpo = log(runif(m, 5, 15)))
    d$y <- rcpb(m, exp(-1 + 0.3 * d$x + d$lexpo + rep(rnorm(length(units), 0, 0.3), each = Tt)), 0.5)
    d$y[rbinom(m, 1, 0.25) == 1] <- 0L
    d
  }
  dh <- mk(1:20); nd <- mk(1001:1003)
  h <- hurdle_cpb(y ~ x, dh, fe = "unit", offset = "lexpo")
  P <- underdisp:::.pmf_and_y_newdata(h, nd, 400)$P
  expect_equal(drop(P %*% (0:400)), as.numeric(predict(h, newdata = nd, type = "response")), tolerance = 1e-8)
  lam <- exp(underdisp:::.fe_ref(h$intensity) + h$int_beta[["x"]] * nd$x + nd$lexpo)
  top <- underdisp:::.cpb_top(lam / (1 - h$intensity$alpha), underdisp:::.cpb_support(h))
  expect_equal(ceiling_exceedance(h, nd)$exceed, sum(nd$y > top))
})
