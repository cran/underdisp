## A unit with no positive count adds nothing to a two-part fit's likelihood, and none of its parameters counts in
## the degrees of freedom (?compare_models, Degrees of freedom): on a panel with such units each class's df equals its
## df on the panel without them, whichever equation carries the unit dummies. The panel's first unit (the reference
## level of the dummies) is one of them.
zero_df_panel <- function(seed = 11, U = 12, TT = 10, zu = c("u01", "u05", "u09")) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(U), function(i) {
    x <- rnorm(TT); a <- rnorm(1, 0.3, 0.4)
    part <- rbinom(TT, 1, plogis(0.4 + 0.5 * x))
    data.frame(unit = sprintf("u%02d", i), t = seq_len(TT), x = x,
               y = part * (1L + rpois(TT, exp(0.5 + a + 0.3 * x))))
  }))
  d$y[d$unit %in% zu] <- 0L
  list(full = d, sub = d[!d$unit %in% zu, ], zu = zu)
}

test_that("the hurdles do not count the participation dummies of units that never participate", {
  p <- zero_df_panel()
  h <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit")
  hs <- hurdle_cpb(y ~ x, p$sub, participation = ~ x, fe = "unit", part_fe = "unit")
  expect_identical(h$df, hs$df)
  expect_setequal(h$part_zero_units, p$zu)
  expect_equal(as.numeric(logLik(h)), as.numeric(logLik(hs)), tolerance = 1e-6)
  expect_identical(attr(logLik(h), "df"), h$df)
  expect_output(print(h), "3 units have no positive count. Their participation dummies are at the limit")
  cm <- suppressMessages(suppressWarnings(compare_models(h = h, z = zi_cpb(y ~ x, p$full, zero = ~ x, fe = "unit",
                                                                         zero_fe = "unit"))))
  expect_equal(cm["h", "df"], h$df)
  ## the floor support counts them, as before 0.2.0
  hf <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit", support = "floor")
  expect_identical(hf$df, length(coef(hf$participation)) + hf$intensity$df)
  expect_length(hf$part_zero_units, 0L)
  ## a hurdle and its count-family twin share the participation equation, so they are charged alike
  hc <- hurdle_count(y ~ x, p$full, family = "poisson", participation = ~ x, fe = "unit", part_fe = "unit",
                     se = "none")
  hcs <- hurdle_count(y ~ x, p$sub, family = "poisson", participation = ~ x, fe = "unit", part_fe = "unit",
                      se = "none")
  expect_identical(hc$df, hcs$df)
  expect_identical(hc$df - hc$intensity$df, h$df - h$intensity$df)
  expect_setequal(hc$part_zero_units, p$zu)
})

test_that("in new data the hurdles give such a unit's rows the participation probability at the limit", {
  p <- zero_df_panel()
  ## a covariate far out: glm's finite dummy for a zero unit plus a large slope term would give such a row a
  ## probability near one
  nd <- p$full[p$full$t == 1, ]; nd$x <- 100; nd$y <- 1L
  zu <- nd$unit %in% p$zu
  h <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit", se = "none")
  ph <- predict(h, newdata = nd, type = "participation")
  expect_identical(ph[zu], rep(0, sum(zu)))
  expect_true(all(ph[!zu] > 0.5))
  ## the held-out scorers read the same probability: a positive count of such a unit is impossible
  P <- underdisp:::.pmf_and_y_newdata(h, nd, kmax = 5L)$P
  expect_identical(P[zu, 1], rep(1, sum(zu)))
  ## fitted values are unchanged (the glm's), and the floor support predicts as before 0.2.0
  expect_identical(h$p_full, as.numeric(fitted(h$participation)))
  hf <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit", se = "none", support = "floor")
  expect_identical(predict(hf, newdata = nd, type = "participation"),
                   as.numeric(predict(hf$participation, newdata = nd, type = "response")))
  ## hurdle_count with intensity dummies too: its intensity never saw such a unit, which enters it at the reference
  hc <- hurdle_count(y ~ x, p$full, family = "poisson", participation = ~ x, fe = "unit", part_fe = "unit", se = "none")
  expect_identical(predict(hc, newdata = nd, type = "participation")[zu], rep(0, sum(zu)))
  expect_identical(underdisp:::.pmf_and_y_newdata(hc, nd, kmax = 5L)$P[zu, 1], rep(1, sum(zu)))
})

test_that("in new data the hurdles give a unit whose every count is positive the participation probability one", {
  p <- zero_df_panel(); d <- p$full; ou <- c("u03", "u07")
  d$y[d$unit %in% ou & d$y == 0L] <- 1L
  ## a covariate far out: glm's finite dummy for such a unit plus a large negative slope term would give its row a
  ## probability near zero
  nd <- d[d$t == 1, ]; nd$x <- -100; nd$y <- 0L
  o <- nd$unit %in% ou
  h <- hurdle_cpb(y ~ x, d, participation = ~ x, fe = "unit", part_fe = "unit", se = "none")
  expect_setequal(h$part_one_units, ou)
  ph <- predict(h, newdata = nd, type = "participation")
  expect_identical(ph[o], rep(1, sum(o)))
  expect_true(all(ph[!o] < 0.5))
  ## the held-out scorers read the same probability: a zero of such a unit is impossible
  expect_identical(underdisp:::.pmf_and_y_newdata(h, nd, kmax = 5L)$P[o, 1], rep(0, sum(o)))
  ## the fit is unchanged: fitted values are the glm's, and the dummies are counted
  expect_identical(h$p_full, as.numeric(fitted(h$participation)))
  expect_identical(h$df, length(coef(h$participation)) - length(h$part_zero_units) + h$intensity$df)
  ## the floor support predicts as before 0.2.0
  hf <- hurdle_cpb(y ~ x, d, participation = ~ x, fe = "unit", part_fe = "unit", se = "none", support = "floor")
  expect_length(hf$part_one_units, 0L)
  expect_identical(predict(hf, newdata = nd, type = "participation"),
                   as.numeric(predict(hf$participation, newdata = nd, type = "response")))
  hc <- hurdle_count(y ~ x, d, family = "poisson", participation = ~ x, part_fe = "unit", se = "none")
  expect_setequal(hc$part_one_units, ou)
  expect_identical(predict(hc, newdata = nd, type = "participation")[o], rep(1, sum(o)))
  expect_identical(underdisp:::.pmf_and_y_newdata(hc, nd, kmax = 5L)$P[o, 1], rep(0, sum(o)))
  g <- hurdle_gec(y ~ x, d, participation = ~ x, part_fe = "unit")
  expect_setequal(g$part_one_units, ou)
  expect_identical(predict(g, newdata = nd, type = "participation")[o], rep(1, sum(o)))
  expect_identical(underdisp:::.pmf_and_y_newdata(g, nd, kmax = 5L)$P[o, 1], rep(0, sum(o)))
})

test_that("hurdle_gec, zi_gec, zi_count and zi_cpb(zero_fe) do not count them either", {
  skip_on_cran()
  p <- zero_df_panel()
  g <- hurdle_gec(y ~ x, p$full, participation = ~ x, part_fe = "unit")
  gs <- hurdle_gec(y ~ x, p$sub, participation = ~ x, part_fe = "unit")
  expect_identical(g$df, gs$df); expect_setequal(g$part_zero_units, p$zu)
  nd <- p$full[p$full$t == 1, ]; nd$x <- 100; nd$y <- 1L; zu <- nd$unit %in% p$zu
  expect_identical(predict(g, newdata = nd, type = "participation")[zu], rep(0, sum(zu)))
  expect_identical(underdisp:::.pmf_and_y_newdata(g, nd, kmax = 5L)$P[zu, 1], rep(1, sum(zu)))
  zg <- suppressWarnings(zi_gec(y ~ x, p$full, zero = ~ 1, zero_fe = "unit"))
  zgs <- suppressWarnings(zi_gec(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit"))
  expect_identical(zg$df, zgs$df); expect_setequal(zg$zero_fe_zero_units, p$zu)
  zc <- suppressWarnings(zi_count(y ~ x, p$full, family = "poisson", zero = ~ 1, fe = "unit", zero_fe = "unit",
                                  se = "none"))
  zcs <- suppressWarnings(zi_count(y ~ x, p$sub, family = "poisson", zero = ~ 1, fe = "unit", zero_fe = "unit",
                                   se = "none"))
  expect_identical(zc$df, zcs$df)
  expect_setequal(zc$fe_zero_units, p$zu); expect_setequal(zc$zero_fe_zero_units, p$zu)
  expect_output(print(zc), "Their count effects and zero-equation dummies are at the limit")
  zc1 <- suppressWarnings(zi_count(y ~ x, p$full, family = "poisson", zero = ~ 1, fe = "unit", se = "none"))
  zc1s <- suppressWarnings(zi_count(y ~ x, p$sub, family = "poisson", zero = ~ 1, fe = "unit", se = "none"))
  expect_identical(zc1$df, zc1s$df)
  ## zi_cpb with zero-equation dummies only: on King's support not counted, on the floor support counted
  z <- suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ 1, zero_fe = "unit"))
  zs <- suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit"))
  expect_identical(z$df, zs$df); expect_setequal(z$zero_fe_zero_units, p$zu)
  zf <- suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ 1, zero_fe = "unit", support = "floor"))
  zfs <- suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit", support = "floor"))
  expect_identical(zf$df, zfs$df + length(p$zu)); expect_length(zf$zero_fe_zero_units, 0L)
  ## with count fixed effects zi_cpb leaves the units out, and the matched hurdle is charged alike
  zfe <- suppressMessages(suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ x, fe = "unit", zero_fe = "unit")))
  zfes <- suppressMessages(suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ x, fe = "unit", zero_fe = "unit")))
  expect_identical(zfe$df, zfes$df)
})
