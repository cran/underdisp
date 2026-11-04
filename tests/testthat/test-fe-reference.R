## Fixed-effects fits: the reference unit of first differences and new rows averages only the units whose effect
## the data identify, and per-row results come back in the data's own row order.

fe_panel <- function(seed = 8) {
  set.seed(seed)
  do.call(rbind, lapply(1:20, function(i) {
    xx <- rnorm(8); data.frame(u = i, t = 1:8, x = xx, y = rcpb(8, exp(1 + rnorm(1, 0, .3) + .3 * xx), .5))
  }))
}

test_that("units with no positive count do not move the reference unit", {
  p <- fe_panel(); p2 <- p; p2$y[p2$u %in% 1:4] <- 0L
  f1 <- cpb_fe(y ~ x, p2, fe = "u"); f2 <- cpb_fe(y ~ x, p2[!p2$u %in% 1:4, ], fe = "u")
  expect_false(any(f1$fe_identified[1:4]))
  expect_equal(f1$fe_ref, f2$fe_ref, tolerance = 1e-4)
  expect_equal(first_difference(f1, "x", -1, 1)$to, first_difference(f2, "x", -1, 1)$to, tolerance = 1e-3)
  g1 <- gec_fe(y ~ x, p2, fe = "u"); g2 <- gec_fe(y ~ x, p2[!p2$u %in% 1:4, ], fe = "u")
  expect_equal(g1$fe_ref, g2$fe_ref, tolerance = 1e-3)
  set.seed(15)
  pc <- do.call(rbind, lapply(1:30, function(i) { xx <- rnorm(12); data.frame(u = i, x = xx, y = rpois(12, exp(.8 + .3 * xx + rnorm(1, 0, .3)))) }))
  pc$y[pc$u %in% 1:4] <- 0L
  r1 <- suppressWarnings(count_reg(y ~ x, pc, family = "poisson", fe = "u"))
  r2 <- count_reg(y ~ x, pc[!pc$u %in% 1:4, ], family = "poisson", fe = "u")
  expect_equal(first_difference(r1, "x", -1, 1)$to, first_difference(r2, "x", -1, 1)$to, tolerance = 1e-4)
})

test_that("under King's zero truncation a unit whose counts are all one stays out of the reference", {
  p <- fe_panel(9); p <- p[p$y > 0, ]; p$y[p$u %in% 1:3] <- 1L
  f <- cpb_fe(y ~ x, p, fe = "u", truncated = TRUE)
  expect_false(any(f$fe_identified[1:3]))
  expect_equal(f$fe_ref, mean(f$fe[f$fe_identified]))
})

test_that("fixed-effects fits return per-row results in the data's row order", {
  p <- fe_panel(); set.seed(2); sh <- p[sample(nrow(p)), ]
  f <- cpb_fe(y ~ x, sh, fe = "u")
  expect_identical(as.numeric(f$Y), as.numeric(sh$y))
  expect_equal(f$loglik, cpb_fe(y ~ x, p, fe = "u")$loglik)
  expect_gt(cor(fitted(f), sh$y), 0.5)
  g <- gec_fe(y ~ x, sh, fe = "u")
  expect_identical(as.numeric(g$Y), as.numeric(sh$y))
  sh$y0 <- ifelse(runif(nrow(sh)) < .3, 0L, sh$y)
  z <- suppressWarnings(zi_cpb(y0 ~ x, sh, zero = ~ 1, fe = "u"))
  expect_identical(as.numeric(z$y), as.numeric(sh$y0))
})

test_that("zi_cpb(fe =) refuses when the unit effects cannot be estimated", {
  set.seed(6); pz <- do.call(rbind, lapply(1:15, function(i) data.frame(u = i, x = rnorm(10))))
  pz$y <- ifelse(runif(150) < .3, 0L, rcpb(150, exp(1 + .3 * pz$x), .5))
  expect_error(zi_cpb(y ~ 1, pz, zero = ~ 1, fe = "u"), "first-stage fixed-effects fit")
})
