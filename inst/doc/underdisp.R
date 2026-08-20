## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)
set.seed(1)

## ----setup--------------------------------------------------------------------
library(underdisp)

## ----sim----------------------------------------------------------------------
n <- 400
x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
y <- rbinom(n, N, 0.5)
d <- data.frame(y = y, x = x)
c(mean = mean(y), var = var(y), ratio = var(y) / mean(y))

## ----screen-------------------------------------------------------------------
ud_screen(y ~ x, data = d, run_cpb = FALSE)

## ----fit----------------------------------------------------------------------
fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
summary(fit)

## ----qoi----------------------------------------------------------------------
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "response")
implied_ceiling(fit, newdata = data.frame(x = 0))

## ----boot---------------------------------------------------------------------
fit_b <- cpb(y ~ x, data = d[d$y > 0, ], se = "bootstrap", B = 99)
summary(fit_b)
irr(fit_b)             # incidence-rate ratios with percentile intervals
alpha_confint(fit_b)   # profile-likelihood interval for alpha

## ----gec----------------------------------------------------------------------
gec(y ~ x, data = d, se = "none")                                  # delta ~ 0.5
gec(y ~ x, data = data.frame(y = rpois(n, exp(1 + 0.4 * x)), x = x),
    se = "none")                                                   # delta ~ 1

## ----fe-----------------------------------------------------------------------
panel <- do.call(rbind, lapply(1:50, function(i) {
  xx <- rnorm(12); NN <- pmax(round(exp(rnorm(1, 0, 0.4) + 0.4 * xx) / 0.5), 1)
  data.frame(unit = i, x = xx, y = rbinom(12, NN, 0.5))
}))
cpb_fe(y ~ x, data = panel, fe = "unit")

## ----family-------------------------------------------------------------------
compare_dispersion(y ~ x, data = d)$table

## ----matched------------------------------------------------------------------
cpb_fit <- cpb(y ~ x, data = d, truncated = FALSE, se = "none")
compare_models(
  CPB          = cpb_fit,
  Poisson      = count_reg(y ~ x, data = d, family = "poisson"),
  NB           = count_reg(y ~ x, data = d, family = "negbin"),
  `COM-Poisson`= count_reg(y ~ x, data = d, family = "compois")
)

## ----pkscreen-----------------------------------------------------------------
data(peacekeeping)
ud_screen(contributions ~ democracy + lgdppc + lpop + milper + factor(iso3),
          data = peacekeeping, run_cpb = FALSE, run_gp = FALSE)

## ----hurdle-------------------------------------------------------------------
z <- rnorm(n)
yh <- rhurdle_cpb(n, lambda = exp(1.2 + 0.3 * x), alpha = 0.5,
                  p = plogis(0.3 + 0.8 * z))
dh <- data.frame(y = yh, x = x, z = z)
h <- hurdle_cpb(y ~ x, data = dh, participation = ~ z)
zi <- zi_cpb(y ~ x, data = dh, zero = ~ z)
compare_models(hurdle = h, mixture = zi)

## ----jackknife----------------------------------------------------------------
short <- do.call(rbind, lapply(1:30, function(i) {
  xx <- rnorm(8); NN <- pmax(round(exp(1.0 + rnorm(1, 0, 0.4) + 0.3 * xx) / 0.5), 1)
  data.frame(unit = i, x = xx, y = rbinom(8, NN, 0.5))
}))
ml <- cpb_fe(y ~ x, data = short, fe = "unit")
jk <- cpb_fe(y ~ x, data = short, fe = "unit", bias_correct = "jackknife")
c(ml = ml$alpha, jackknife = jk$alpha)   # truth is 0.5; ML is biased downward

## ----dharma, eval = requireNamespace("DHARMa", quietly = TRUE)----------------
sims <- simulate(h, nsim = 100, seed = 1)
res <- DHARMa::createDHARMa(simulatedResponse = as.matrix(sims),
                            observedResponse  = dh$y,
                            fittedPredictedResponse = fitted(h),
                            integerResponse = TRUE)
plot(res)

