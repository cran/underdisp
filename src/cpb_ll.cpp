// cpb_ll.cpp -- the pooled (zero-truncated) CPB negative log-likelihood and the
// per-observation log pmf. Exact model (matches R/cpb.R): support {0,...,floor(n_i)},
// n_i = lambda_i/(1-alpha). The pmf is normalized by the mode-outward recursion
// of cpb_pmf.h, and the pure-R likelihood (.cpb_pmf_core) is its independent check.
#include <Rcpp.h>
#include <vector>
#include "cpb_pmf.h"
using namespace Rcpp;

// [[Rcpp::export]]
double cpb_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                   int max_support, bool truncated) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));       // logit^-1
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return 1e10;
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double eta = offset[i];
    for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double contrib = cpb_logpmf(Y[i], std::exp(eta), alpha, la, l1a, max_support, truncated);
    if (contrib <= -1e299 || !R_finite(contrib)) return 1e10;   // above its ceiling, or the alpha ~ 1 guard
    ll += contrib;
  }
  return -ll;
}

// Weighted CPB negative log-likelihood, -sum_i w_i log P(Y_i), for the ZI-CPB
// M-step (the count component is fit with weights 1 - responsibility).
// [[Rcpp::export]]
double cpb_wnll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y,
                    NumericVector w, NumericVector offset, int max_support, bool truncated) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return 1e10;
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double eta = offset[i]; for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double contrib = cpb_logpmf(Y[i], std::exp(eta), alpha, la, l1a, max_support, truncated);
    if (contrib <= -1e299 || !R_finite(contrib)) return 1e10;
    ll += w[i] * contrib;
  }
  return -ll;
}

// Per-observation log P(Y = y_i) (column 0) and P(Y = 0) (column 1) for the CPB,
// for the ZI-CPB E-step and observed-data log-likelihood. params = (beta, logit alpha).
// [[Rcpp::export]]
NumericMatrix cpb_lp0_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                          int max_support, bool truncated) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; i++) {
    double eta = offset[i]; for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double p0, lp = cpb_logpmf_p0(Y[i], std::exp(eta), alpha, la, l1a, max_support, truncated, p0);   // -1e300 where infeasible
    out(i, 0) = lp; out(i, 1) = p0;
  }
  return out;
}
