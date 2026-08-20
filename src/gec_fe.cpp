// gec_fe.cpp -- King GEC / Katz-family regression with high-dimensional unit
// fixed effects, concentrated (profiled) out of the likelihood exactly as in
// cpb_fe.cpp. Each unit's intercept is a 1-D golden-section inner maximization;
// the outer optimizer sees only (beta, log delta). Data SORTED BY UNIT.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static double gec_logpmf_fe(int y, double mu, double delta, int max_support) {
  double gamma = (delta - 1.0) / delta, omega = mu / delta;
  std::vector<double> lw; lw.push_back(0.0);
  double mx = 0.0; int k = 0;
  while (k <= max_support) {
    double num = omega + gamma * (double)k;
    if (num <= 0.0) break;
    double next = lw[k] + std::log(num) - std::log((double)(k + 1));
    lw.push_back(next); if (next > mx) mx = next; k++;
    if (gamma >= 0.0 && next - mx < -36.0) break;
  }
  int K = (int)lw.size() - 1;
  if (y > K) return -1e300;
  double s = 0.0; for (int j = 0; j <= K; j++) s += std::exp(lw[j] - mx);
  return lw[y] - (mx + std::log(s));
}

static inline double unit_ll(double a, const IntegerVector& Y, const std::vector<double>& o,
                             int lo, int hi, double delta, int max_support) {
  double ll = 0.0;
  for (int t = lo; t < hi; t++) {
    double lp = gec_logpmf_fe(Y[t], std::exp(a + o[t]), delta, max_support);
    if (lp <= -1e299) return -1e300;
    ll += lp;
  }
  return ll;
}

static double maximize_unit(const IntegerVector& Y, const std::vector<double>& o, int lo, int hi,
                            double delta, int max_support, double& a_star, int inner_it) {
  double alo = -1e300; int ymax = 0;
  for (int t = lo; t < hi; t++) if (Y[t] > 0) {
    if (Y[t] > ymax) ymax = Y[t];
    if (delta < 1.0) {                              // finite support: mu/(1-delta) >= Y
      double b = std::log((1.0 - delta) * (double)Y[t]) - o[t];
      if (b > alo) alo = b;
    }
  }
  if (ymax == 0) { a_star = -30.0; return 0.0; }
  double L, H;
  if (delta < 1.0) { L = alo + 1e-6; H = alo + 25.0; } else { L = -20.0; H = 20.0; }
  const double gr = 0.6180339887498949;
  double c = H - gr * (H - L), d = L + gr * (H - L);
  double fc = unit_ll(c, Y, o, lo, hi, delta, max_support);
  double fd = unit_ll(d, Y, o, lo, hi, delta, max_support);
  for (int it = 0; it < inner_it; it++) {
    if (fc < fd) { L = c; c = d; fc = fd; d = L + gr * (H - L);
                   fd = unit_ll(d, Y, o, lo, hi, delta, max_support); }
    else        { H = d; d = c; fd = fc; c = H - gr * (H - L);
                   fc = unit_ll(c, Y, o, lo, hi, delta, max_support); }
  }
  a_star = 0.5 * (L + H);
  return unit_ll(a_star, Y, o, lo, hi, delta, max_support);
}

static void compute_offsets(NumericVector params, NumericMatrix X, NumericVector uoff, std::vector<double>& o) {
  int n = X.nrow(), p = X.ncol();
  for (int i = 0; i < n; i++) { double e = uoff[i]; for (int j = 0; j < p; j++) e += X(i, j) * params[j]; o[i] = e; }
}

// [[Rcpp::export]]
double gec_fe_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                      IntegerVector ustart, int n_units, int max_support, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  double delta = std::exp(params[p]);
  if (delta <= 1e-8 || delta >= 1e8) return 1e10;
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  double total = 0.0, astar;
  for (int u = 0; u < n_units; u++) {
    double ll = maximize_unit(Y, o, ustart[u], ustart[u + 1], delta, max_support, astar, inner_it);
    if (ll <= -1e299) return 1e10;
    total += ll;
  }
  return -total;
}

// [[Rcpp::export]]
NumericVector gec_fe_intercepts_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                                    IntegerVector ustart, int n_units, int max_support, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  double delta = std::exp(params[p]);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector fe(n_units); double astar;
  for (int u = 0; u < n_units; u++) {
    maximize_unit(Y, o, ustart[u], ustart[u + 1], delta, max_support, astar, inner_it);
    fe[u] = astar;
  }
  return fe;
}
