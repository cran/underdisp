## The one rule for the uncertainty of a quantity of interest (.qoi_ci): joint bootstrap draws (coefficients,
## dispersion, reference unit), recentered for fixed-effects fits, else the delta method over the full covariance;
## `at =` profiles; the shared bootstrap layout; predictions for seen units; and the refusals around them.

sim_panel <- function(seed, N = 30, TT = 8) {          # 30 units: the cluster bootstrap's lower guard
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) { x <- rnorm(TT)
    data.frame(unit = i, t = seq_len(TT), x = x, y = rcpb(TT, exp(0.8 + rnorm(1, 0, 0.5) + 0.4 * x), 0.5)) }))
}

test_that("the bootstrap draws carry the dispersion and, for fixed effects, the reference unit", {
  skip_on_cran()
  set.seed(1); n <- 250; x <- rnorm(n); d <- data.frame(y = rcpb(n, exp(1.3 + 0.4 * x), 0.5), x = x)
  f <- suppressMessages(cpb(y ~ x, d, truncated = FALSE, se = "bootstrap", B = 25))
  expect_identical(colnames(f$boot), c("(Intercept)", "x", "alpha"))
  expect_identical(attr(f$boot, "B"), 25L); expect_identical(attr(f$boot, "cluster"), "obs")
  dp <- sim_panel(2)
  fe <- suppressMessages(cpb_fe(y ~ x, dp, fe = "unit", se = "bootstrap", B = 20))
  expect_identical(colnames(fe$boot), c("x", "alpha", "fe_ref"))
  expect_identical(attr(fe$boot, "cluster"), "unit")
  expect_equal(unname(sqrt(diag(vcov(fe)))), unname(fe$se.beta), tolerance = 1e-12)
  g <- suppressMessages(gec_fe(y ~ x, dp, fe = "unit", se = "bootstrap", B = 15))
  expect_identical(colnames(g$boot), c("x", "delta", "fe_ref"))
})

test_that("a fixed-effects first difference recenters the draws on the estimate and moves alpha with them", {
  skip_on_cran()
  dp <- sim_panel(3)
  fe <- suppressMessages(cpb_fe(y ~ x, dp, fe = "unit", se = "bootstrap", B = 30))
  fd <- first_difference(fe, "x", -1, 1)
  expect_identical(fd$method, "bootstrap (stored)")
  expect_true(fd$lower < fd$diff && fd$diff < fd$upper)
  ## by hand: the draws recentered, the CPB mean at each
  D <- fe$boot[complete.cases(fe$boot), ]; th <- c(fe$coefficients, alpha = fe$alpha, fe_ref = fe$fe_ref)
  D <- sweep(D, 2, colMeans(D) - th)
  q <- apply(D, 1, function(r) diff(underdisp:::.cpb_mean(exp(r[["fe_ref"]] + c(-1, 1) * r[["x"]]), r[["alpha"]],
                                                         FALSE, "ceiling")))
  expect_equal(c(fd$lower, fd$upper), unname(quantile(q, c(0.025, 0.975))), tolerance = 1e-8)
})

test_that("the delta method of the count families includes the dispersion parameter", {
  set.seed(4); n <- 300; x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(y = rgammacount(n, exp(1 + 0.4 * x), alpha = 2), x = x, z = z)
  m <- count_reg(y ~ x + z, d, family = "gammacount")
  fd <- first_difference(m, "x", 0, 1)
  fam <- underdisp:::.count_fam("gammacount"); th <- c(m$coefficients, log(m$theta))
  f <- function(p) { xr <- c(1, 0, mean(z)); mu <- exp(c(sum(xr * p[1:3]), sum(replace(xr, 2, 1) * p[1:3])))
                     diff(fam$meanfun(mu, exp(p[4]))) }
  g <- numDeriv::grad(f, th); se <- sqrt(drop(t(g) %*% m$vcov_full %*% g))
  expect_equal(fd$upper - fd$diff, qnorm(0.975) * se, tolerance = 1e-4)
  expect_identical(fd$method, "delta")
  ## the zero and marginal rows of zi_count now have intervals
  set.seed(5); yz <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
  zc <- zi_count(y ~ x, data.frame(y = yz, x = x, z = z), family = "poisson", zero = ~ z)
  fz <- first_difference(zc, "z", -1, 1)
  expect_true(all(is.finite(c(fz$lower, fz$upper)))); expect_identical(fz$method, rep("delta", 3L))
})

test_that("`at` fixes other columns and is checked", {
  set.seed(6); n <- 200; x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(y = rpois(n, exp(0.8 + 0.3 * x + 0.2 * z)), x = x, z = z)
  m <- count_reg(y ~ x + z, d, family = "poisson")
  fd <- first_difference(m, "x", 0, 1, at = c(z = 2))
  expect_equal(c(fd$from, fd$to), predict(m, newdata = data.frame(x = c(0, 1), z = 2)), tolerance = 1e-10)
  expect_error(first_difference(m, "x", 0, 1, at = c(w = 1)), "not in the model")
  expect_error(first_difference(m, "x", 0, 1, at = c(x = 1)), "'variable' moves")
  expect_error(first_difference(m, "x", 0, 1, at = 1), "named numeric")
})

test_that("a hurdle's participation equation is clustered with `cluster`, and its first differences use it", {
  skip_if_not_installed("sandwich")
  set.seed(7); nn <- 500; cl <- rep(1:50, each = 10); x <- rnorm(nn); z <- rnorm(nn) + rep(rnorm(50), each = 10)
  d <- data.frame(y = rhurdle_cpb(nn, exp(1.2 + 0.4 * x), 0.5, plogis(-0.2 + 0.8 * z)), x = x, z = z, cl = cl)
  h <- hurdle_cpb(y ~ x, d, participation = ~ z, cluster = "cl")
  expect_identical(h$part_vcov_type, "cluster")
  expect_equal(unname(h$part_vcov), unname(sandwich::vcovCL(h$participation, cluster = d$cl, type = "HC0")),   # glm's bread is from its last IRLS weights
               tolerance = 1e-3)
  fd <- first_difference(h, "z", -1, 1)
  expect_identical(fd$method[1], "delta"); expect_identical(fd$method[3], "none")   # no intensity inference
  hc <- hurdle_count(y ~ x, d, family = "poisson", participation = ~ z, se = "robust")
  expect_equal(unname(hc$part_vcov), unname(sandwich::vcovHC(hc$participation, type = "HC0")), tolerance = 1e-3)
})

test_that("a separated participation equation is flagged, and its ratios and standard errors are withheld", {
  set.seed(8); n <- 300; x <- rnorm(n); s <- rbinom(n, 1, 0.3)
  p <- ifelse(s == 1, 1, plogis(0.2 + 0.5 * x))                     # every s = 1 unit participates
  d <- data.frame(y = ifelse(rbinom(n, 1, p) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L), x = x, s = s)
  w <- character(0)
  h <- withCallingHandlers(hurdle_count(y ~ x, d, family = "poisson", participation = ~ x + s),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("participation equation is separated", w)))
  expect_identical(h$part_separation$diverging, "s")
  ir <- suppressMessages(irr(h))
  expect_true(is.na(ir$estimate[ir$term == "s" & ir$equation == "participation"]))
  expect_true(is.na(summary(h)$blocks[[1]]["s", "Std. Error"]))
  cm <- suppressWarnings(compare_models(a = h, b = hurdle_count(y ~ x, d, family = "poisson", participation = ~ x)))
  expect_true(cm["a", "separated"]); expect_false(cm["b", "separated"])
})

test_that("zi_test refuses the cases its boundary reference does not cover", {
  skip_on_cran()
  set.seed(9); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, lambda = exp(1.2 + 0.3 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z)); d <- data.frame(y = y, x = x, z = z)
  nest <- suppressWarnings(cpb(y ~ x, d, truncated = FALSE))   # at the Poisson limit on zero-inflated data
  expect_error(zi_test(zi_cpb(y ~ x, d, zero = ~ z), nest), "covariate-dependent")
  expect_s3_class(zi_test(zi_cpb(y ~ x, d, zero = ~ 1), nest), "zi_test")
})

test_that("predict() on fixed-effects fits uses a seen unit's own effect, and gives no interval there", {
  dp <- sim_panel(10)
  fe <- cpb_fe(y ~ x, dp, fe = "unit")
  expect_equal(predict(fe, newdata = dp), unname(fitted(fe)), tolerance = 1e-10)
  nd <- dp[1:3, ]; nd$unit <- 999
  expect_equal(predict(fe, newdata = nd, type = "rate"), exp(fe$fe_ref + nd$x * fe$coefficients[["x"]]), tolerance = 1e-10)
  expect_equal(predict(fe, newdata = nd[, c("x", "y")], type = "rate"), exp(fe$fe_ref + nd$x * fe$coefficients[["x"]]),
               tolerance = 1e-10)
  expect_message(r <- predict(fe, newdata = dp[1:2, ], se.fit = TRUE), "own fixed effect")
  expect_true(all(is.na(r$se.fit))); expect_identical(r$method, rep("none", 2L))
  expect_error(predict(fe, se.fit = TRUE), "need 'newdata'")
})

test_that("predict() intervals and ame() follow the first difference's rule on every class", {
  skip_on_cran()
  set.seed(21); n <- 300; x <- rnorm(n); z <- rnorm(n)
  d0 <- data.frame(y = rcpb(n, exp(1.2 + 0.5 * x), 0.5), x = x, z = z, u = factor(sample(1:30, n, TRUE)))
  ## zero-inflated data for the mixtures, so that their zero equations move the mean
  dz <- transform(d0, y = rzicpb(n, exp(1.2 + 0.5 * x), 0.5, plogis(-1 + 0.8 * z)))
  ## overdispersed data for the hurdle negative binomial: on d0 its theta runs to the Poisson limit (log theta
  ## near 18), where theta no longer moves the mean. Drawn on its own seed with the stream restored, so every
  ## other fit's data and bootstrap draws are as before
  dn <- local({ old <- .Random.seed; on.exit(assign(".Random.seed", old, envir = globalenv()))
                set.seed(23); transform(d0, y = stats::rnbinom(n, mu = exp(1.2 + 0.5 * x), size = 2)) })
  q <- function(e) suppressWarnings(suppressMessages(e))
  fits <- list(cpb = q(cpb(y ~ x + z, d0, truncated = FALSE, se = "bootstrap", B = 20)),
               cpb_fe = q(cpb_fe(y ~ x + z, d0, fe = "u", se = "bootstrap", B = 15)),
               gec = q(gec(y ~ x + z, d0, se = "bootstrap", B = 20)),
               gec_fe = q(gec_fe(y ~ x + z, d0, fe = "u", se = "bootstrap", B = 15)),
               hurdle_cpb = q(hurdle_cpb(y ~ x + z, d0, participation = ~ x + z, se = "bootstrap", B = 20)),
               hurdle_gec = q(hurdle_gec(y ~ x + z, d0, participation = ~ x + z, se = "bootstrap", B = 20)),
               zi_cpb = q(zi_cpb(y ~ x + z, dz, zero = ~ x + z, se = "bootstrap", B = 15)),
               zi_gec = q(zi_gec(y ~ x + z, dz, zero = ~ x + z, se = "bootstrap", B = 15)),
               count_reg = q(count_reg(y ~ x + z, d0, family = "compois")),   # its mean moves with nu
               hurdle_count = q(hurdle_count(y ~ x + z, dn, family = "negbin", participation = ~ x + z)),
               zi_count = q(zi_count(y ~ x + z, dz, family = "poisson", zero = ~ x + z)))
  nd <- d0[1:5, c("x", "z")]; nd0 <- data.frame(x = 0, z = mean(z))
  qoi_par <- function(f) eval(quote(.qoi_par(f)), list(f = f), asNamespace("underdisp"))   # an internal generic
  for (nm in names(fits)) {
    f <- fits[[nm]]; p <- qoi_par(f); pt <- predict(f, nd)
    ## the parameters enter the prediction, each of them, and at the estimate it is the fit's own
    expect_equal(predict(p$set(p$theta), nd), pt, tolerance = 1e-12, label = paste(nm, "identity"))
    for (k in names(p$theta)) { th <- p$theta; th[k] <- th[k] + 0.05
      expect_gt(max(abs(predict(p$set(th), nd) - pt)), 1e-10, label = paste(nm, k)) }
    pr <- q(predict(f, nd, se.fit = TRUE, interval = "confidence"))
    expect_equal(pr$fit[, "fit"], pt, label = nm)
    expect_true(all(pr$fit[, "lwr"] < pr$fit[, "upr"] & pr$se.fit > 0), label = paste(nm, "interval"))
    ## ame() over one row at the reference profile is the first difference, interval included
    fd <- q(first_difference(f, "x", 0, 1)); r <- nrow(fd)
    am <- q(ame(f, "x", nd0, from = 0, to = 1))
    expect_equal(c(am$estimate, am$lower, am$upper), c(fd$diff[r], fd$lower[r], fd$upper[r]), tolerance = 1e-6,
                 label = paste(nm, "ame = fd"))
    expect_identical(am$method, fd$method[r])
  }
  ## the derivative, by central differences over the rows
  f <- fits$cpb; h <- 1e-4 * diff(range(d0$x))
  hand <- mean((predict(f, transform(d0, x = x + h / 2)) - predict(f, transform(d0, x = x - h / 2))) / h)
  a1 <- ame(f, "x", d0); expect_equal(a1$estimate, hand, tolerance = 1e-10); expect_identical(a1$contrast, "dY/dx")
  ## fixed effects: rows at their own unit effects have no interval; at the reference unit they do
  expect_message(a2 <- ame(fits$cpb_fe, "x", d0), "own fixed effect"); expect_true(is.na(a2$lower))
  a3 <- ame(fits$cpb_fe, "x", d0[, c("x", "z")]); expect_identical(a3$method, "bootstrap (stored)")
  ## a factor contrast through an interaction, and the refusals
  d0$g <- factor(ifelse(d0$z > 0, "b", "a"))
  fg <- q(count_reg(y ~ x * g, d0, family = "negbin"))
  ag <- ame(fg, "g", d0, from = "a", to = "b")
  hand <- mean(predict(fg, transform(d0, g = factor("b", levels = c("a", "b")))) -
               predict(fg, transform(d0, g = factor("a", levels = c("a", "b")))))
  expect_equal(ag$estimate, hand, tolerance = 1e-12); expect_identical(ag$method, "delta")
  expect_error(ame(fg, "g", d0), "numeric variable")
  expect_error(ame(fg, "g", d0, from = "a", to = "c"), "not a level")
  expect_error(ame(fg, "u", d0), "not a variable of the model")
  expect_error(ame(fits$cpb_fe, "u", d0), "unit column")
  expect_error(ame(fg, "x", d0, from = 0), "both 'from' and 'to'")
})

test_that("zi_cpb(fe) leaves units with no positive count out of the estimation", {
  skip_on_cran()
  set.seed(7); U <- 30; TT <- 6
  dp <- do.call(rbind, lapply(1:U, function(i) { x <- rnorm(TT)
    data.frame(unit = sprintf("u%02d", i), x = x, y = rzicpb(TT, exp(0.6 + rnorm(1, 0, .5) + 0.4 * x), 0.5, 0.3)) }))
  zu <- sprintf("u%02d", c(3, 11, 17, 25)); dp$y[dp$unit %in% zu] <- 0L
  dp <- dp[sample(nrow(dp)), ]; sub <- dp[!dp$unit %in% zu, ]; r <- dp$unit %in% zu
  for (m in c("ml", "em")) {
    zf <- if (m == "em") "unit"
    expect_message(f <- zi_cpb(y ~ x, dp, zero = ~ 1, fe = "unit", zero_fe = zf, method = m), "4 units (24 rows)", fixed = TRUE)
    g <- suppressMessages(zi_cpb(y ~ x, sub, zero = ~ 1, fe = "unit", zero_fe = zf, method = m))
    ## the limit: the same fit as on the panel without those units, with every row kept and aligned
    expect_equal(f$coefficients, g$coefficients, tolerance = 1e-12); expect_equal(f$alpha, g$alpha, tolerance = 1e-12)
    expect_equal(f$zero_coef, g$zero_coef, tolerance = 1e-12); expect_equal(f$loglik, g$loglik, tolerance = 1e-12)
    expect_identical(f$df, g$df); expect_identical(f$n, nrow(dp)); expect_identical(f$y, as.integer(dp$y))
    expect_setequal(f$zero_units, zu); expect_true(all(fitted(f)[r] < 1e-8)); expect_false(anyNA(fitted(f)))
    expect_true(all(simulate(f, 2, seed = 1)[r, ] == 0))
    expect_equal(predict(f, data.frame(unit = "u03", x = 0)), 0, tolerance = 1e-8)
    expect_true(is.finite(score(f)[["logscore"]]))
  }
  expect_true(all(is.na(f$pi_full[r])))                         # zero_fe: their inflation dummy is not identified
  expect_true(is.na(predict(f, data.frame(unit = "u03", x = 0), type = "zero")))
  ## the floor support keeps the fill of versions before 0.2.0
  fl <- zi_cpb(y ~ x, dp, zero = ~ 1, fe = "unit", support = "floor")
  expect_null(fl$n_zero_units)
  expect_false(isTRUE(all.equal(fl$loglik, zi_cpb(y ~ x, sub, zero = ~ 1, fe = "unit", support = "floor")$loglik)))
})

test_that("the alpha-limit columns of implied_ceiling(), and alpha's interval without a bootstrap", {
  set.seed(11); n <- 200; x <- rnorm(n); d <- data.frame(y = rcpb(n, exp(1.4 + 0.3 * x), 0.5), x = x)
  f <- cpb(y ~ x, d, truncated = FALSE)
  ic <- implied_ceiling(f, newdata = data.frame(x = 0))
  expect_true(all(c("ceiling_at_alpha_lower", "ceiling_at_alpha_upper") %in% names(ic)))
  expect_false(any(c("lower", "upper") %in% names(ic)))
  ci <- confint(f, "alpha"); expect_true(ci[1] < f$alpha && f$alpha < ci[2])
  expect_error(confint(f), "refit with se")
})

test_that("Pearson residuals, the dispersion's direction in glance(), and tidy() intervals", {
  skip_if_not_installed("broom")
  set.seed(12); n <- 250; x <- rnorm(n); d <- data.frame(y = rcpb(n, exp(1.4 + 0.3 * x), 0.5), x = x)
  g <- gec(y ~ x, d); m <- count_reg(y ~ x, d, family = "gammacount")
  nb <- count_reg(y ~ x, data.frame(y = rnbinom(n, mu = exp(1 + 0.3 * x), size = 2), x = x), family = "negbin")
  for (f in list(g, m)) {
    r <- residuals(f, type = "pearson"); expect_true(all(is.finite(r)))
    expect_equal(sign(r), sign(residuals(f)))
  }
  expect_identical(broom::glance(cpb(y ~ x, d, truncated = FALSE))$dispersion_direction, "underdispersed")
  expect_identical(broom::glance(m)$dispersion_direction, "underdispersed")
  expect_identical(broom::glance(m)$dispersion_parameter, "alpha")
  expect_identical(broom::glance(nb)$dispersion_direction, "overdispersed")
  td <- broom::tidy(m, conf.int = TRUE)
  expect_equal(unname(as.matrix(td[, c("conf.low", "conf.high")])), unname(confint(m)[td$term, ]), tolerance = 1e-12)
})

test_that("small fixes: a vector cluster follows sorted rows, `|` inside I(), a family object, two-way fe", {
  skip_on_cran()
  dp <- sim_panel(13); set.seed(1); sh <- dp[sample(nrow(dp)), ]
  sh$grp <- sh$unit + 100                                       # 30 clusters, given in shuffled order
  set.seed(3); a <- suppressMessages(cpb_fe(y ~ x, sh, fe = "unit", se = "bootstrap", B = 10, cluster = "grp"))
  set.seed(3); b <- suppressMessages(cpb_fe(y ~ x, sh, fe = "unit", se = "bootstrap", B = 10, cluster = sh$grp))
  expect_equal(a$boot, b$boot, ignore_attr = TRUE)
  set.seed(14); d <- data.frame(y = rpois(100, 2), a = rbinom(100, 1, 0.5) == 1, b = rbinom(100, 1, 0.5) == 1)
  expect_s3_class(count_reg(y ~ I(a | b), d, family = "poisson"), "count_reg")
  expect_error(count_reg(y ~ a | b, d), "two-part formula")
  expect_equal(count_reg(y ~ a, d, family = poisson)$loglik, count_reg(y ~ a, d, family = "poisson")$loglik)
  expect_error(cpb_fe(y ~ x, dp, fe = c("unit", "t")), "two-way effects")
  ## a covariate called alpha does not collide with the dispersion parameter
  set.seed(15); al <- rnorm(200); dd <- data.frame(y = rcpb(200, exp(1.3 + 0.4 * al), 0.5), alpha = al)
  fa <- suppressMessages(cpb(y ~ alpha, dd, truncated = FALSE, se = "bootstrap", B = 20))
  fd <- first_difference(fa, "alpha", -1, 1)
  expect_identical(fd$method, "bootstrap (stored)"); expect_true(fd$lower < fd$diff && fd$diff < fd$upper)
})

test_that("delta-method intervals do not depend on a covariate's units", {
  set.seed(1); n <- 500; pop <- runif(n, 1e5, 9e5); z <- rnorm(n)
  y <- rpois(n, exp(-0.5 + 1.2 * pop / 1e6 + 0.3 * z))
  fm <- count_reg(y ~ pop + z, data.frame(y, z, pop = pop / 1e6), family = "poisson", se = "robust")
  fp <- count_reg(y ~ pop + z, data.frame(y, z, pop = pop), family = "poisson", se = "robust")
  a <- first_difference(fm, "z", 0, 1); b <- first_difference(fp, "z", 0, 1)
  expect_equal(c(b$lower, b$upper), c(a$lower, a$upper), tolerance = 1e-6)
  expect_equal(predict(fp, data.frame(pop = 5e5, z = 0), se.fit = TRUE)$se.fit,
               predict(fm, data.frame(pop = 0.5, z = 0), se.fit = TRUE)$se.fit, tolerance = 1e-6)
})

test_that("ame() weights its average, and refuses the floor support's slope", {
  set.seed(2); tab <- data.frame(x = rep(c(-1, 0, 1), 4), y = rpois(12, 3)); tab$w <- sample(1:9, 12, TRUE)
  fw <- count_reg(y ~ x, tab, family = "poisson", weights = "w")
  ex <- tab[rep(seq_len(nrow(tab)), tab$w), ]
  expect_equal(ame(fw, "x", tab, weights = "w")$estimate, ame(count_reg(y ~ x, ex, family = "poisson"), "x", ex)$estimate,
               tolerance = 1e-6)
  expect_message(ame(fw, "x", tab), "frequency weights")
  expect_error(ame(fw, "x", tab, weights = -tab$w), "non-negative")
  set.seed(3); xx <- rnorm(200)
  ff <- suppressWarnings(cpb(y ~ x, data.frame(y = rcpb(200, exp(1 + 0.3 * xx), 0.5, support = "floor"), x = xx),
                             truncated = FALSE, support = "floor"))
  expect_error(ame(ff, "x", data.frame(x = xx)), "jumps at integer ceilings")
  expect_s3_class(ame(ff, "x", data.frame(x = xx), from = 0, to = 1), "ud_ame")
})

test_that("ame(B =) refits on resampled units and averages each replicate over its own units", {
  skip_on_cran()
  set.seed(3)
  pd <- do.call(rbind, lapply(1:30, function(i) {
    x <- rnorm(6); data.frame(unit = i, x = x, y = rcpb(6, exp(rnorm(1, 0, 0.5) + 0.4 * x), 0.5))
  }))
  ff <- cpb_fe(y ~ x, data = pd, fe = "unit")
  expect_message(a0 <- ame(ff, "x", pd), "ame\\(B = \\)")
  expect_identical(a0$method, "none")
  set.seed(1); a1 <- suppressMessages(ame(ff, "x", pd, B = 20))
  expect_identical(a1$method, "bootstrap (refit, B=20)")
  expect_equal(a1$estimate, a0$estimate)
  expect_true(a1$lower < a1$estimate && a1$estimate < a1$upper)
  expect_error(ame(ff, "x", pd[, c("unit", "x")], B = 5), "lacks 'y'")
  expect_error(ame(ff, "x", pd, B = 1), "at least 2")
  ## a weights vector given to the fit follows the resampled rows
  set.seed(4); x <- rnorm(150); w <- rep(1:3, 50); dd <- data.frame(y = rpois(150, exp(0.5 + 0.3 * x)), x = x)
  cr <- count_reg(y ~ x, dd, family = "poisson", se = "robust", weights = w)
  set.seed(2); a2 <- suppressMessages(ame(cr, "x", dd, weights = w, B = 30))
  expect_identical(a2$method, "bootstrap (refit, B=30)")
  expect_true(a2$lower < a2$estimate && a2$estimate < a2$upper)
})

test_that("ame(B =) follows the fit's own jackknife outcome, and ignores a factor's unused levels", {
  skip_on_cran()
  set.seed(202)
  d <- do.call(rbind, lapply(1:20, function(i) {             # a dispersion regime change: the gate refuses
    x <- rnorm(12); lam <- exp(1.0 + rnorm(1, 0, 0.4) + 0.3 * x); a <- ifelse(1:12 <= 6, 0.3, 0.7)
    data.frame(unit = i, tt = 1:12, x = x, y = vapply(1:12, function(t) rcpb(1, lam[t], a[t]), integer(1)))
  }))
  jk <- suppressWarnings(cpb_fe(y ~ x, data = d, fe = "unit", bias_correct = "jackknife", time = "tt"))
  expect_identical(jk$bias_correct, "none")
  raw <- cpb_fe(y ~ x, data = d, fe = "unit")
  set.seed(11); a1 <- suppressMessages(ame(jk, "x", d, from = -1, to = 1, B = 6))
  set.seed(11); a2 <- suppressMessages(ame(raw, "x", d, from = -1, to = 1, B = 6))
  expect_equal(c(a1$lower, a1$upper), c(a2$lower, a2$upper))
  expect_error(ame(jk, "x", d[, c("unit", "x", "y")], B = 2), "lacks 'tt'")
  d2 <- d; d2$unit <- factor(d2$unit, levels = 1:30)        # ten unused levels
  set.seed(2); b1 <- suppressMessages(ame(raw, "x", d2, B = 6))
  d3 <- d2; d3$unit <- droplevels(d3$unit)
  set.seed(2); b2 <- suppressMessages(ame(raw, "x", d3, B = 6))
  expect_equal(c(b1$lower, b1$upper), c(b2$lower, b2$upper))
})

test_that("a recentred draw whose dispersion leaves its range is dropped", {
  D <- cbind(b = c(0.1, 0.2, 0.3, 0.4, 0.5), .disp = c(0.80, 0.84, 0.90, 0.70, 0.76))
  th <- c(b = 0.3, .disp = 0.95)                 # the draws' mean is 0.80: recentring adds 0.15
  seen <- numeric(0)
  q <- function(t) { seen <<- c(seen, t[[".disp"]]); t[[".disp"]] }
  ci <- underdisp:::.qoi_ci(q, th, draws = D, recenter = TRUE, bounds = list(.disp = c(0, 1)))
  kept <- D[, ".disp"] - (mean(D[, ".disp"]) - 0.95); kept <- kept[kept <= 1]
  expect_true(all(seen <= 1)); expect_length(kept, 4L)
  expect_equal(ci$upper, unname(stats::quantile(kept, 0.975)))
  expect_gt(underdisp:::.qoi_ci(q, th, draws = D, recenter = TRUE)$upper, 1)     # unbounded, it enters
  skip_on_cran()
  fe <- cpb_fe(y ~ x, sim_panel(41), fe = "unit")
  expect_identical(underdisp:::.qoi_par(fe)$bounds, list(.disp = c(0, 1)))
})

test_that("draws saved before version 0.2.0 are read in their own layout", {
  skip_on_cran()
  set.seed(18); d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rcpb(200, exp(1 + 0.3 * d$x), 0.5)
  ## rbind() of the replicates' named vectors: the coefficients' names, an empty one for the dispersion, and none
  ## of the attributes .ud_boot() sets
  old_layout <- function(D, keep, nm) { D <- unclass(D)[, keep, drop = FALSE]; dimnames(D) <- list(NULL, nm); D }
  g <- suppressMessages(gec(y ~ x, d, se = "bootstrap", B = 5))
  og <- g; og$boot <- old_layout(g$boot, colnames(g$boot), c(names(g$coefficients), ""))
  expect_null(attr(og$boot, "B"))
  expect_equal(confint(og), confint(g))
  d$y[rbinom(200, 1, 0.3) == 1] <- 0L
  z <- suppressMessages(zi_cpb(y ~ x, d, zero = ~ z, se = "bootstrap", B = 5))
  cz <- c(paste0("count:", names(z$coefficients)), paste0("zero:", names(z$zero_coef)))
  oz <- z; oz$boot <- old_layout(z$boot, cz, c(names(z$coefficients), names(z$zero_coef)))   # no alpha column
  expect_equal(vcov(oz), vcov(z))
})

test_that("a bootstrap with fewer than two converged replicates gives NA intervals", {
  skip_on_cran()
  set.seed(18); d <- data.frame(x = rnorm(200)); d$y <- rcpb(200, exp(1 + 0.3 * d$x), 0.5)
  g <- suppressMessages(gec(y ~ x, d, se = "bootstrap", B = 5)); g$boot[-1, ] <- NA
  expect_warning(cg <- confint(g), "fewer than two")
  expect_identical(rownames(cg), c("(Intercept)", "x", "delta")); expect_true(all(is.na(cg)))
  f <- suppressMessages(cpb(y ~ x, d, truncated = FALSE, se = "bootstrap", B = 5)); f$boot[-1, ] <- NA
  expect_warning(cf <- confint(f), "fewer than two")
  expect_true(all(is.na(cf[c("(Intercept)", "x"), ]))); expect_true(all(is.finite(cf["alpha", ])))
})
