#' underdisp: Diagnostics and Models for Underdispersed Count Data
#'
#' Detect and model underdispersion (conditional variance below the conditional
#' mean) in count data. The package provides a screening diagnostic
#' ([ud_screen()]), the continuous parameter binomial regression model
#' ([cpb()]) with a zero-truncated variant, validated bootstrap and
#' profile-likelihood inference, a family-comparison helper
#' ([compare_dispersion()]), and King-style quantities of interest
#' ([predict.cpb()], [implied_ceiling()], [first_difference()]).
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
