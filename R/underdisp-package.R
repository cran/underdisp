#' underdisp: Diagnostics and Models for Underdispersed Count Data
#'
#' Detect and model underdispersion (conditional variance below the conditional
#' mean) in count data.
#'
#' @section Screening and diagnostics:
#' [ud_screen()] (the zero-truncated-Poisson at-risk screen with calibrated or
#' parametric-bootstrap thresholds), [dispersion_test()] (regression-adjusted
#' tests of equidispersion), [dispersion_profile()] (the conditional
#' variance-to-mean curve against each family's implied curve),
#' [compare_dispersion()], [zi_test()], [rootogram()], [pit_hist()].
#'
#' @section Estimators:
#' The hard-ceiling family [cpb()] / [cpb_fe()] and the free-dispersion family
#' [gec()] / [gec_fe()], with hurdle and zero-inflated forms
#' ([hurdle_cpb()], [zi_cpb()], [hurdle_gec()], [zi_gec()]); the matched count
#' families through one interface, [count_reg()], [hurdle_count()],
#' [zi_count()] (Poisson, negative binomial, COM-Poisson in the rate and the
#' mean parameterization, generalized Poisson, gamma-count, double Poisson),
#' all with offsets, frequency weights, fixed effects, and analytic, robust,
#' cluster, or bootstrap inference.
#'
#' @section Comparison and quantities of interest:
#' [compare_models()], [score()], [cv_score()]; [predict()],
#' [implied_ceiling()], [irr()], [first_difference()] (with the
#' extensive/intensive decomposition for the two-part models), [confint()] for
#' every class, [simulate()] for DHARMa diagnostics, and broom / texreg /
#' modelsummary support.
#'
#' @keywords internal
#' @useDynLib underdisp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.frame model.matrix model.response delete.response terms reformulate
#' @importFrom stats glm glm.fit poisson binomial quasibinomial optim uniroot
#' @importFrom stats nobs predict confint fitted residuals coef logLik AIC vcov simulate
#' @importFrom stats dpois dnbinom plogis qlogis pnorm qnorm pchisq qchisq ppois printCoefmat
#' @importFrom stats complete.cases na.omit cov var sd quantile setNames rnorm rbinom update
#' @importFrom graphics plot rect lines points abline barplot
"_PACKAGE"
