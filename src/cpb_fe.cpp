// cpb_fe.cpp -- CPB regression with high-dimensional unit fixed effects via a
// CONCENTRATED (profiled) likelihood. Each unit's intercept is solved by a 1-D inner
// maximization, so the outer optimizer sees only (beta, alpha). No unit dummies.
// Data must be SORTED BY UNIT; ustart[u]..ustart[u+1]-1 are unit u's rows.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// log P_CPB(Y=y | lambda, alpha), using a precomputed log-factorial table lf.
static inline double cpb_logpmf(int y, double lam, double alpha, double la, double l1a,
                                const std::vector<double>& lf, int max_support, bool truncated) {
  double ni = lam / (1.0 - alpha);
  int Ki = (int)std::floor(ni);
  if (y > Ki || Ki > max_support) return -1e300;
  double lgni1 = R::lgammafn(ni + 1.0);
  double g = lgni1, mx = -1e300;
  // first pass: max
  std::vector<double> lw(Ki + 1);
  for (int k = 0; k <= Ki; k++) {
    if (k > 0) g -= std::log(ni - k + 1.0);
    lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la;
    if (lw[k] > mx) mx = lw[k];
  }
  double s = 0.0; for (int k = 0; k <= Ki; k++) s += std::exp(lw[k] - mx);
  double logD = mx + std::log(s);
  double lp = lw[y] - logD;
  if (truncated) {
    double e = std::exp(lw[0] - logD);
    if (e >= 1.0 - 1e-15) return -1e300;
    lp -= std::log1p(-e);
  }
  return lp;
}

// unit log-likelihood at intercept a over contiguous rows [lo,hi).
static inline double unit_ll(double a, const IntegerVector& Y, const std::vector<double>& o,
                             int lo, int hi, double alpha, double la, double l1a,
                             const std::vector<double>& lf, int max_support, bool truncated) {
  double ll = 0.0;
  for (int t = lo; t < hi; t++) {
    double lp = cpb_logpmf(Y[t], std::exp(a + o[t]), alpha, la, l1a, lf, max_support, truncated);
    if (lp <= -1e299) return -1e300;
    ll += lp;
  }
  return ll;
}

// maximize unit_ll over the intercept a (golden section on a feasible bracket).
static double maximize_unit(const IntegerVector& Y, const std::vector<double>& o, int lo, int hi,
                            double alpha, double la, double l1a, const std::vector<double>& lf,
                            int max_support, bool truncated, double& a_star, int inner_it) {
  double alo = -1e300; int ymax = 0;
  for (int t = lo; t < hi; t++)
    if (Y[t] > 0) {
      double b = std::log(Y[t] * (1.0 - alpha)) - o[t];
      if (b > alo) alo = b;
      if (Y[t] > ymax) ymax = Y[t];
    }
  if (ymax == 0) { a_star = -30.0; return 0.0; }
  double L = alo + 1e-6, H = alo + 25.0;
  const double gr = 0.6180339887498949;
  double c = H - gr * (H - L), d = L + gr * (H - L);
  double fc = unit_ll(c, Y, o, lo, hi, alpha, la, l1a, lf, max_support, truncated);
  double fd = unit_ll(d, Y, o, lo, hi, alpha, la, l1a, lf, max_support, truncated);
  for (int it = 0; it < inner_it; it++) {
    if (fc < fd) { L = c; c = d; fc = fd; d = L + gr * (H - L);
                   fd = unit_ll(d, Y, o, lo, hi, alpha, la, l1a, lf, max_support, truncated); }
    else        { H = d; d = c; fd = fc; c = H - gr * (H - L);
                   fc = unit_ll(c, Y, o, lo, hi, alpha, la, l1a, lf, max_support, truncated); }
  }
  a_star = 0.5 * (L + H);
  return unit_ll(a_star, Y, o, lo, hi, alpha, la, l1a, lf, max_support, truncated);
}

static void compute_offsets(NumericVector params, NumericMatrix X, NumericVector uoff, std::vector<double>& o) {
  int n = X.nrow(), p = X.ncol();
  for (int i = 0; i < n; i++) { double e = uoff[i]; for (int j = 0; j < p; j++) e += X(i, j) * params[j]; o[i] = e; }
}

// [[Rcpp::export]]
double cpb_fe_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                      IntegerVector ustart, int n_units, int max_support,
                      bool truncated, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return 1e10;
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  double total = 0.0, astar;
  for (int u = 0; u < n_units; u++) {
    double ll = maximize_unit(Y, o, ustart[u], ustart[u+1], alpha, la, l1a, lf, max_support, truncated, astar, inner_it);
    if (ll <= -1e299) return 1e10;
    total += ll;
  }
  return -total;
}

// [[Rcpp::export]]
NumericVector cpb_fe_intercepts_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                                    IntegerVector ustart, int n_units, int max_support,
                                    bool truncated, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector fe(n_units); double astar;
  for (int u = 0; u < n_units; u++) {
    maximize_unit(Y, o, ustart[u], ustart[u+1], alpha, la, l1a, lf, max_support, truncated, astar, inner_it);
    fe[u] = astar;
  }
  return fe;
}
