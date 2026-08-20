## Regression tests for the QoI / predict / parity defects found in the v0.1.0
## adversarial re-review (predict recycling, dropped offsets, truncated prob,
## zero_fe predict, GEC/FE QoI parity). Each guards a specific fixed bug.

test_that("predict.hurdle_cpb marginal E(Y) is full-length and not recycled", {
  set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.6 * z)) == 1, 0L, 1L + rpois(n, exp(0.6 + 0.4 * x)))
  hc <- hurdle_cpb(y ~ x, data.frame(y = y, x = x, z = z), participation = ~ z, se = "none")
  pm <- predict(hc, type = "response")
  expect_length(pm, n)
  expect_equal(unname(pm), unname(hc$p_full * hc$lambda_full), tolerance = 1e-9)  # not the recycled short vector
})

test_that("predict.cpb applies the stored offset", {
  set.seed(2); n <- 300; x <- rnorm(n); e <- log(runif(n, 1, 4))
  d <- data.frame(y = rpois(n, exp(e + 0.5 * x)), x = x, le = e)
  cf <- cpb(y ~ x, d[d$y > 0, ], truncated = FALSE, se = "none", offset = "le")
  expect_equal(unname(predict(cf, type = "response")), unname(fitted(cf)), tolerance = 1e-8)
  expect_equal(unname(predict(cf, type = "link")), unname(cf$linear.predictors), tolerance = 1e-8)
})

test_that("predict(type='prob') conditions on Y>0 for truncated fits", {
  set.seed(3); n <- 400
  ct <- count_reg(y ~ 1, data.frame(y = 1L + rpois(n, 2)), truncated = TRUE, se = "none")
  mu <- exp(unname(ct$coefficients[1]))
  expect_equal(predict(ct, data.frame(x = 0), type = "prob", at = 1), underdisp:::dztpois(1, mu), tolerance = 1e-6)
  gt <- gec(y ~ x, data.frame(y = 1L + rpois(n, exp(0.5 + 0.3 * rnorm(n))), x = rnorm(n)),
            truncated = TRUE, se = "none")
  expect_equal(predict(gt, data.frame(x = 0), type = "prob", at = 0), 0)
})

test_that("zero_fe predict works on a single-level newdata (xlev stored)", {
  skip_on_cran()
  set.seed(4); n <- 400; x <- rnorm(n); z <- rnorm(n); g <- factor(sample(5, n, replace = TRUE))
  y <- ifelse(rbinom(n, 1, plogis(-0.4 + 0.7 * z)) == 1, 0L, rzicpb(n, exp(1.1 + 0.4 * x), 0.5, 0))
  d <- data.frame(y = y, x = x, z = z, g = g)
  nd <- data.frame(x = 0, z = 0, g = factor("1", levels = levels(g)))
  expect_silent(predict(zi_cpb(y ~ x, d, zero = ~ z, zero_fe = "g", se = "none"), nd, type = "zero"))
  expect_silent(predict(zi_gec(y ~ x, d, zero = ~ z, zero_fe = "g", se = "none"), nd, type = "zero"))
})

test_that("GEC irr/first_difference expose bootstrap CIs and the FE/two-part QoIs exist", {
  skip_on_cran()
  set.seed(5); n <- 300; x <- rnorm(n); z <- rnorm(n); g <- factor(sample(6, n, replace = TRUE))
  d <- data.frame(y = rpois(n, exp(0.6 + 0.5 * x)), x = x)
  gb <- gec(y ~ x, d, se = "bootstrap", B = 25)
  expect_true(all(c("lower", "upper") %in% colnames(irr(gb))))
  expect_true(all(c("lower", "upper") %in% colnames(first_difference(gb, "x", 0, 1))))
  N <- pmax(round(exp(1.4 + 0.4 * x) / 0.5), 1); yb <- rbinom(n, N, 0.5); df <- data.frame(y = yb, x = x, g = g)
  expect_s3_class(first_difference(gec_fe(y ~ x, df, fe = "g", se = "none"), "x", 0, 1), "data.frame")
  expect_s3_class(first_difference(cpb_fe(y ~ x, df, fe = "g", se = "none"), "x", 0, 1), "data.frame")
  expect_no_error(irr(cpb_fe(y ~ x, df, fe = "g", se = "none")))
  yh <- ifelse(rbinom(n, 1, plogis(0.3 + 0.6 * z)) == 1, 0L, 1L + rpois(n, exp(0.6 + 0.4 * x)))
  dh <- data.frame(y = yh, x = x, z = z)
  expect_equal(nrow(first_difference(hurdle_gec(y ~ x, dh, participation = ~ z, se = "none"), "x", 0, 1)), 3L)
  dz <- data.frame(y = ifelse(rbinom(n, 1, plogis(-0.4 + 0.7 * z)) == 1, 0L, rpois(n, exp(1.1 + 0.4 * x))),
                   x = x, z = z)
  expect_equal(nrow(first_difference(zi_gec(y ~ x, dz, zero = ~ z, se = "none"), "x", 0, 1)), 3L)
})

test_that("summary() gives a coefficient display for the FE/two-part/GEC classes", {
  set.seed(6); n <- 300; x <- rnorm(n); g <- factor(sample(6, n, replace = TRUE))
  N <- pmax(round(exp(1.4 + 0.4 * x) / 0.5), 1)
  df <- data.frame(y = rbinom(n, N, 0.5), x = x, g = g)
  expect_s3_class(summary(cpb_fe(y ~ x, df, fe = "g", se = "none")), "summary.underdisp")
  expect_s3_class(summary(gec(y ~ x, df, se = "none")), "summary.underdisp")
})

test_that("coef() of a ZI fit returns both blocks (count + zero), like pscl", {
  set.seed(7); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(-0.4 + 0.7 * z)) == 1, 0L, rpois(n, exp(1.1 + 0.4 * x)))
  cz <- coef(zi_count(y ~ x, data.frame(y = y, x = x, z = z), zero = ~ z))
  expect_named(cz, c("count", "zero"))
})
