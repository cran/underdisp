// cpb_ll.cpp -- fast C++ evaluation of the (zero-truncated) CPB negative log-likelihood.
// Exact model (matches cpb.R / cpb_zt.R): support {0,...,floor(n_i)}, n_i = lambda_i/(1-alpha).
// Speed comes from a gamma recurrence: lgamma(n_i - k + 1) is computed once (k=0) and then
// updated as g_k = g_{k-1} - log(n_i - k + 1), so only ONE lgamma call per observation.
// Validated against the pure-R likelihood to machine precision (see 10_cpb_cpp_validate.R).
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
double cpb_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                   int max_support, bool truncated) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));       // logit^-1
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return 1e10;
  double la = std::log(alpha), l1a = std::log1p(-alpha);

  // precompute log(k!) up to max_support
  std::vector<double> lf(max_support + 2);
  lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);

  double ll = 0.0;
  std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    double eta = offset[i];
    for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double lam = std::exp(eta);
    double ni  = lam / (1.0 - alpha);
    int Ki = (int)std::floor(ni);
    if (Y[i] > Ki) return 1e10;            // observation exceeds implied ceiling
    if (Ki > max_support) return 1e10;     // alpha ~ 1 guard

    double lgni1 = R::lgammafn(ni + 1.0);
    double g = lgni1;                       // g_k = lgamma(ni - k + 1); g_0 = lgamma(ni+1)
    double mx = -1e300;
    lw.assign(Ki + 1, 0.0);
    for (int k = 0; k <= Ki; k++) {
      if (k > 0) g -= std::log(ni - k + 1.0);
      lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la;
      if (lw[k] > mx) mx = lw[k];
    }
    double s = 0.0;
    for (int k = 0; k <= Ki; k++) s += std::exp(lw[k] - mx);
    double logD = mx + std::log(s);

    double contrib = lw[Y[i]] - logD;                 // log P(Y_i = y_i)
    if (truncated) {                                  // condition on Y >= 1
      double e = std::exp(lw[0] - logD);              // P(Y_i = 0)
      if (e >= 1.0 - 1e-15) return 1e10;
      contrib -= std::log1p(-e);
    }
    if (!R_finite(contrib)) return 1e10;
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
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  double ll = 0.0; std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    double eta = offset[i]; for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double lam = std::exp(eta), ni = lam / (1.0 - alpha);
    int Ki = (int)std::floor(ni);
    if (Y[i] > Ki || Ki > max_support) return 1e10;
    double lgni1 = R::lgammafn(ni + 1.0), g = lgni1, mx = -1e300;
    lw.assign(Ki + 1, 0.0);
    for (int k = 0; k <= Ki; k++) { if (k > 0) g -= std::log(ni - k + 1.0);
      lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la; if (lw[k] > mx) mx = lw[k]; }
    double s = 0.0; for (int k = 0; k <= Ki; k++) s += std::exp(lw[k] - mx);
    double logD = mx + std::log(s);
    double contrib = lw[Y[i]] - logD;
    if (truncated) { double e = std::exp(lw[0] - logD); if (e >= 1.0 - 1e-15) return 1e10;
      contrib -= std::log1p(-e); }
    if (!R_finite(contrib)) return 1e10;
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
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  NumericMatrix out(n, 2); std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    double eta = offset[i]; for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double lam = std::exp(eta), ni = lam / (1.0 - alpha);
    int Ki = (int)std::floor(ni);
    if (Ki > max_support || Y[i] > Ki) { out(i, 0) = -1e300; out(i, 1) = 0.0; continue; }
    double lgni1 = R::lgammafn(ni + 1.0), g = lgni1, mx = -1e300;
    lw.assign(Ki + 1, 0.0);
    for (int k = 0; k <= Ki; k++) { if (k > 0) g -= std::log(ni - k + 1.0);
      lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la; if (lw[k] > mx) mx = lw[k]; }
    double s = 0.0; for (int k = 0; k <= Ki; k++) s += std::exp(lw[k] - mx);
    double logD = mx + std::log(s);
    double lp = lw[Y[i]] - logD, p0 = std::exp(lw[0] - logD);
    if (truncated) { if (p0 >= 1.0 - 1e-15) lp = -1e300; else lp -= std::log1p(-p0); }
    out(i, 0) = lp; out(i, 1) = p0;
  }
  return out;
}
