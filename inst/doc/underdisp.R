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
first_difference(fit, "x", from = -1, to = 1)

## ----dtest--------------------------------------------------------------------
dispersion_test(fit, B = 0)   # asymptotic value; the default is the bootstrap p-value
dispersion_test(count_reg(y ~ x, data = d, family = "poisson"),
                method = "auxiliary", alternative = "under")

## ----boot---------------------------------------------------------------------
fit_b <- cpb(y ~ x, data = d[d$y > 0, ], se = "bootstrap", B = 49)   # a small B keeps the vignette quick
summary(fit_b)
irr(fit_b)             # rate ratios with percentile intervals
confint(fit_b)         # coefficients (percentile) and alpha (first-order profile likelihood)

## ----gec----------------------------------------------------------------------
gec(y ~ x, data = d, se = "none")                                  # delta ~ 0.5
gec(y ~ x, data = data.frame(y = rpois(n, exp(1 + 0.4 * x)), x = x),
    se = "none")                                                   # delta ~ 1

## ----fe-----------------------------------------------------------------------
panel <- do.call(rbind, lapply(1:30, function(i) {
  xx <- rnorm(12); NN <- pmax(round(exp(rnorm(1, 0, 0.4) + 0.4 * xx) / 0.5), 1)
  data.frame(unit = i, x = xx, y = rbinom(12, NN, 0.5))
}))
fe_fit <- cpb_fe(y ~ x, data = panel, fe = "unit")
fe_fit
dispersion_test(fe_fit, B = 0)   # the statistic against a Poisson with the same unit effects

## ----fe-boot, eval = FALSE----------------------------------------------------
# dispersion_test(fe_fit, cores = 2)   # parametric-bootstrap p-value, 199 replicates

## ----family-------------------------------------------------------------------
compare_dispersion(y ~ x, data = d)$table

## ----wider--------------------------------------------------------------------
dw <- d[1:150, ]                         # part of the sample keeps the seven fits quick
fits <- list(
  CPB          = cpb(y ~ x, data = dw, truncated = FALSE, se = "none"),
  Poisson      = count_reg(y ~ x, data = dw, family = "poisson"),
  NB           = count_reg(y ~ x, data = dw, family = "negbin"),
  `COM-Poisson`= count_reg(y ~ x, data = dw, family = "mpcmp"),
  GenPoisson   = count_reg(y ~ x, data = dw, family = "genpois"),
  GammaCount   = count_reg(y ~ x, data = dw, family = "gammacount"),
  DoublePois   = count_reg(y ~ x, data = dw, family = "doublepois")
)
do.call(compare_models, fits)

## ----profile------------------------------------------------------------------
dispersion_profile(CPB = fits$CPB, NB = fits$NB, GammaCount = fits$GammaCount,
                   GenPoisson = fits$GenPoisson)

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

## ----weights------------------------------------------------------------------
w <- sample(1:3, n, TRUE)
c(weighted = count_reg(y ~ x, data = d, family = "gammacount", weights = w)$loglik,
  expanded = count_reg(y ~ x, data = d[rep(seq_len(n), w), ], family = "gammacount")$loglik)

## ----jackknife----------------------------------------------------------------
short <- do.call(rbind, lapply(1:16, function(i) {
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

