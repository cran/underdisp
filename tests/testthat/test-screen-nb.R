## ud_screen(): the NB-vs-Poisson p-value follows the boundary rule in bootstrap mode

test_that("the NB p-value is asymptotic in the fast mode and bootstrapped in bootstrap mode", {
  set.seed(21); x <- rnorm(150)
  d <- data.frame(y = rpois(150, exp(0.8 + 0.3 * x)), x = x)
  s0 <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE))
  expect_identical(s0$p_nb_method, "asymptotic")
  expect_equal(s0$p_nb, 0.5 * pchisq(max(s0$lr_nb, 0), 1, lower.tail = FALSE))
  expect_output(print(s0), "asymptotic, conservative")
  set.seed(5)
  s1 <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE,
                                   ztp_threshold = "bootstrap", ztp_boot_B = 19))
  after1 <- runif(1)
  set.seed(5)
  s2 <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE,
                                   ztp_threshold = "bootstrap", ztp_boot_B = 19, nb_boot = FALSE))
  after2 <- runif(1)
  expect_identical(s1$p_nb_method, "parametric bootstrap")
  expect_identical(s2$p_nb_method, "asymptotic")
  expect_true(s1$p_nb >= 1 / 20 && s1$p_nb <= 1)
  expect_output(print(s1), "by parametric bootstrap")
  ## the NB bootstrap leaves the random-number stream exactly where it found it
  expect_identical(s1$ztp_threshold, s2$ztp_threshold)
  expect_identical(s1$pearson_ztp, s2$pearson_ztp)
  expect_identical(after1, after2)
})

test_that("the NB bootstrap respects the comparator gates", {
  set.seed(22); x <- rnorm(150)
  d <- data.frame(y = rpois(150, exp(0.8 + 0.3 * x)), x = x)
  s <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE,
                                  ztp_threshold = "bootstrap", ztp_boot_B = 19, comp_max_par = 1))
  expect_identical(s$p_nb_method, "asymptotic")
})

test_that("an NB fit at the Poisson boundary reports LR = 0 and no overdispersion reading", {
  set.seed(23); x <- rnorm(200)
  d <- data.frame(y = rcpb(200, exp(1.2 + 0.3 * x), 0.4, truncated = FALSE), x = x)
  s <- suppressWarnings(ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE))
  expect_identical(s$lr_nb, 0)
  out <- capture.output(print(s))
  expect_true(any(grepl("NB vs Poisson LR = 0, at the Poisson boundary", out, fixed = TRUE)))
  expect_false(any(grepl("sig => overdispersion", out, fixed = TRUE)))
})
