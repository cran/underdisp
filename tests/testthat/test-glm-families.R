## Matched Poisson/NB families: validate against the reference implementations
## (glm, glm.nb, VGAM, sandwich, pscl) and check the score()/compare_models()
## integration. Reference packages beyond Imports are skipped when unavailable.

test_that("count_reg Poisson matches glm", {
  set.seed(42); n <- 500; x <- rnorm(n)
  y <- rpois(n, exp(0.7 + 0.6 * x)); d <- data.frame(y = y, x = x)
  m <- count_reg(y ~ x, d, family = "poisson", se = "analytic")
  g <- glm(y ~ x, d, family = poisson())
  expect_equal(unname(coef(m)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(m$loglik, as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(m$se.beta), unname(sqrt(diag(vcov(g)))), tolerance = 1e-4)
  expect_equal(m$df, length(coef(g)))
})

test_that("count_reg negbin matches glm.nb", {
  set.seed(7); n <- 600; x <- rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(0.6 + 0.5 * x), theta = 2); d <- data.frame(y = y, x = x)
  m <- count_reg(y ~ x, d, family = "negbin")
  g <- suppressWarnings(MASS::glm.nb(y ~ x, d))
  expect_equal(unname(coef(m)), unname(coef(g)), tolerance = 1e-3)
  expect_equal(m$theta, g$theta, tolerance = 1e-2)
  expect_equal(m$loglik, as.numeric(logLik(g)), tolerance = 1e-4)
})

test_that("zero-truncated Poisson matches VGAM::pospoisson", {
  skip_if_not_installed("VGAM")
  set.seed(9); n <- 500; x <- rnorm(n)
  y <- rpois(n, exp(0.9 + 0.4 * x)); d <- data.frame(y = y, x = x); d <- d[d$y > 0, ]
  m <- count_reg(y ~ x, d, family = "poisson", truncated = TRUE)
  v <- suppressWarnings(VGAM::vglm(y ~ x, VGAM::pospoisson, data = d))
  expect_equal(unname(coef(m)), unname(coef(v)), tolerance = 1e-3)
  expect_equal(m$loglik, as.numeric(VGAM::logLik(v)), tolerance = 1e-3)
})

test_that("fixed effects match glm with factor dummies and count the FE df", {
  set.seed(11); n <- 700; x <- rnorm(n); g <- factor(sample(15, n, replace = TRUE))
  y <- rpois(n, exp(0.5 + 0.5 * x)); d <- data.frame(y = y, x = x, g = g)
  m <- count_reg(y ~ x, d, family = "poisson", fe = "g")
  gg <- glm(y ~ x + factor(g), d, family = poisson())
  expect_equal(unname(coef(m)["x"]), unname(coef(gg)["x"]), tolerance = 1e-5)
  expect_equal(m$loglik, as.numeric(logLik(gg)), tolerance = 1e-5)
  expect_equal(m$df, length(coef(gg)))
})

test_that("robust and cluster SE match sandwich", {
  skip_if_not_installed("sandwich")
  set.seed(13); n <- 600; x <- rnorm(n); cl <- sample(30, n, replace = TRUE)
  y <- rpois(n, exp(0.6 + 0.5 * x)); d <- data.frame(y = y, x = x)
  g <- glm(y ~ x, d, family = poisson())
  mr <- count_reg(y ~ x, d, family = "poisson", se = "robust")
  expect_equal(unname(mr$se.beta["x"]),
               unname(sqrt(diag(sandwich::vcovHC(g, "HC0")))["x"]), tolerance = 1e-3)
  mc <- count_reg(y ~ x, d, family = "poisson", se = "cluster", cluster = cl)
  expect_equal(unname(mc$se.beta["x"]),
               unname(sqrt(diag(sandwich::vcovCL(g, cluster = cl, cadjust = FALSE)))["x"]),
               tolerance = 1e-2)
})

test_that("zi_count and hurdle_count match pscl", {
  skip_if_not_installed("pscl")
  set.seed(11); n <- 1200; x <- rnorm(n); z <- rnorm(n)
  pistar <- plogis(-0.5 + 0.9 * z); mu <- exp(0.6 + 0.5 * x)
  y <- ifelse(rbinom(n, 1, pistar) == 1, 0L, rpois(n, mu)); d <- data.frame(y = y, x = x, z = z)
  mzi <- zi_count(y ~ x, d, family = "poisson", zero = ~ z)
  pz <- pscl::zeroinfl(y ~ x | z, data = d, dist = "poisson")
  expect_equal(unname(coef(mzi)$count["x"]), unname(pz$coefficients$count["x"]), tolerance = 1e-3)
  expect_equal(mzi$loglik, as.numeric(logLik(pz)), tolerance = 1e-3)
  mh <- hurdle_count(y ~ x, d, family = "poisson", participation = ~ z)
  ph <- pscl::hurdle(y ~ x | z, data = d, dist = "poisson", zero.dist = "binomial")
  expect_equal(mh$loglik, as.numeric(logLik(ph)), tolerance = 1e-3)
})

test_that("score() and compare_models() accept the new families and rank sensibly", {
  set.seed(3); n <- 900; x <- rnorm(n); z <- rnorm(n)
  pistar <- plogis(-0.4 + 0.8 * z); mu <- exp(0.7 + 0.4 * x)
  y <- ifelse(rbinom(n, 1, pistar) == 1, 0L, rpois(n, mu)); d <- data.frame(y = y, x = x, z = z)
  mp <- count_reg(y ~ x, d, family = "poisson")
  mzi <- zi_count(y ~ x, d, family = "poisson", zero = ~ z)
  mh <- hurdle_count(y ~ x, d, family = "poisson", participation = ~ z)
  expect_named(score(mp), c("logscore", "rps"))
  expect_true(all(is.finite(score(mzi))))
  tab <- compare_models(pois = mp, zip = mzi, hurdle = mh)
  expect_true(all(c("df", "logLik", "AIC", "logscore") %in% names(tab)))
  # on zero-inflated data the mixture/hurdle beat the plain Poisson on AIC
  expect_lt(tab["zip", "AIC"], tab["pois", "AIC"])
})

test_that("predict returns response and link scales", {
  set.seed(5); n <- 300; x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.5 * x)); d <- data.frame(y = y, x = x)
  m <- count_reg(y ~ x, d, family = "poisson")
  nd <- data.frame(x = c(-1, 0, 1))
  expect_equal(predict(m, nd, type = "response"), exp(predict(m, nd, type = "link")))
  expect_length(predict(m, nd), 3L)
})

test_that("zero-truncated negbin does not strand below its Poisson limit on underdispersed data", {
  ## Regression: on tight (underdispersed) positives the ZT-NB likelihood is
  ## monotone in theta toward the Poisson limit; a single theta = 1 start left
  ## BFGS stranded ~900 log-likelihood points below the ZT-Poisson it nests.
  set.seed(11); n <- 600; x <- rnorm(n)
  y <- rcpb(n, exp(1 + 0.4 * x), 0.55, truncated = TRUE)
  d <- data.frame(y = y, x = x)
  pz <- count_reg(y ~ x, d, family = "poisson", truncated = TRUE, se = "none")
  nb <- count_reg(y ~ x, d, family = "negbin",  truncated = TRUE, se = "none")
  expect_gte(nb$loglik, pz$loglik - 0.1)       # NB nests Poisson at theta -> Inf
})

test_that("truncated compois never fits worse than the ZT-Poisson it nests", {
  skip_on_cran()
  ## Regression: on the peacekeeping dynamic-model data the ZT-COM-Poisson BFGS
  ## stalled exactly at its nu = 1 start (the ZT-Poisson), ~78 log-likelihood
  ## points below its profile optimum near nu = 1.5. The second dispersion start
  ## must lift the fit off the nest point when the data are genuinely tight.
  set.seed(12); n <- 800; x <- rnorm(n); w <- rpois(n, 3) + 1   # lag-like covariate
  y <- rcpb(n, exp(0.3 + 0.25 * x + 0.15 * w), 0.6, truncated = TRUE)
  d <- data.frame(y = y, x = x, w = w)
  pz <- count_reg(y ~ x + w, d, family = "poisson", truncated = TRUE, se = "none")
  co <- count_reg(y ~ x + w, d, family = "compois", truncated = TRUE, se = "none")
  expect_gte(co$loglik, pz$loglik - 1e-6)      # COMP nests Poisson at nu = 1
  expect_gt(co$theta, 1.1)                     # and must detect the tightness
})
