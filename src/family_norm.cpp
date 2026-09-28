// family_norm.cpp -- normalizing constants and exact moments of the double
// Poisson (Efron 1986) and the generalized Poisson (Consul & Jain 1973) with
// constant dispersion, one pass per observation. The terms and the support rule
// mirror R/count_families.R (.dp_lterm, .gp_lterm, .grow_support) exactly; the
// R functions delegate here so that a fit does not loop over observations in R.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// log(k!) table, extended on demand
static const std::vector<double>& lgam_tab(int K) {
  static std::vector<double> lg(1, 0.0);
  if ((int)lg.size() <= K) {
    size_t old = lg.size(); lg.resize((size_t)K + 1);
    for (size_t k = old; k <= (size_t)K; k++) lg[k] = lg[k - 1] + std::log((double)k);
  }
  return lg;
}

// ---- double Poisson: log term of the pmf (R: .dp_lterm) ---------------------
static inline double dp_lterm(int k, double mu, double theta, double lmu, double lth, double lgk1) {
  double base = 0.5 * lth - theta * mu;
  if (k == 0) return base;
  double lk = std::log((double)k);
  return base + (-(double)k + k * lk - lgk1) + theta * k * (1.0 + lmu - lk);
}

// support rule (R: .grow_support with K0 = ceil(mu + 12 sqrt(mu/theta + 1) + 30))
static int dp_kmax(double mu, double theta, double lmu, double lth, int cap) {
  double K0 = std::ceil(mu + 12.0 * std::sqrt(mu / theta + 1.0) + 30.0);
  if (!R_finite(K0)) return cap;
  int K = std::max(50, (int)std::min(K0, (double)cap));
  int step = std::max(16, (int)std::min(std::ceil(K * 0.25), 10000.0));
  const std::vector<double>& lg = lgam_tab(K + step + 1);
  double mx = -1e300;
  for (int k = 0; k <= K; k++) { double t = dp_lterm(k, mu, theta, lmu, lth, lg[k]); if (R_finite(t) && t > mx) mx = t; }
  while (K < cap) {
    const std::vector<double>& lg2 = lgam_tab(K + step + 1);
    double tk = dp_lterm(K + step, mu, theta, lmu, lth, lg2[K + step]);
    if (!R_finite(tk) || tk < mx - 40.0) break;
    K += step;
  }
  return std::min(K + step, cap);
}

// n x 3: log Z, mean, variance of the double Poisson at each mu
// [[Rcpp::export]]
NumericMatrix dp_norm_cpp(NumericVector mu, double theta, int cap) {
  int n = mu.size(); NumericMatrix out(n, 3);
  double lth = std::log(theta);
  for (int i = 0; i < n; i++) {
    double m = mu[i];
    if (!(m > 0.0) || !R_finite(m) || !(theta > 0.0)) { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; out(i, 2) = NA_REAL; continue; }
    double lmu = std::log(m);
    int K = dp_kmax(m, theta, lmu, lth, cap);
    const std::vector<double>& lg = lgam_tab(K);
    double mx = -1e300;
    std::vector<double> lt(K + 1);
    for (int k = 0; k <= K; k++) { lt[k] = dp_lterm(k, m, theta, lmu, lth, lg[k]); if (lt[k] > mx) mx = lt[k]; }
    double s = 0.0, s1 = 0.0, s2 = 0.0;
    for (int k = 0; k <= K; k++) { double w = std::exp(lt[k] - mx); s += w; s1 += k * w; s2 += (double)k * k * w; }
    double e1 = s1 / s;
    out(i, 0) = mx + std::log(s); out(i, 1) = e1; out(i, 2) = s2 / s - e1 * e1;
  }
  return out;
}

// ---- generalized Poisson: log term (R: .gp_lterm) ---------------------------
static inline double gp_lterm(int k, double mu, double lambda, double lgk1) {
  double th = mu * (1.0 - lambda), a = th + lambda * k;
  if (a <= 0.0) return R_NegInf;
  return std::log(th) + (k - 1.0) * std::log(std::max(a, 1e-300)) - th - lambda * k - lgk1;
}

static int gp_kmax(double mu, double lambda, int cap) {
  // the terms are unimodal in k for every lambda, so the tail rule applies; for
  // lambda < 0 the finite support end caps the search (the terms beyond are -Inf)
  if (lambda < 0.0) {
    double kend = std::ceil(mu * (1.0 - lambda) / (-lambda)) - 1.0;
    if (kend < (double)cap) cap = std::max(0, (int)kend);
  }
  double K0 = std::ceil(mu + 12.0 * std::sqrt(mu / ((1.0 - lambda) * (1.0 - lambda)) + 1.0) + 30.0);
  if (!R_finite(K0)) return cap;
  int K = std::min(std::max(50, (int)std::min(K0, (double)cap)), cap);
  int step = std::max(16, (int)std::min(std::ceil(K * 0.25), 10000.0));
  const std::vector<double>& lg = lgam_tab(K + step + 1);
  double mx = -1e300;
  for (int k = 0; k <= K; k++) { double t = gp_lterm(k, mu, lambda, lg[k]); if (R_finite(t) && t > mx) mx = t; }
  while (K < cap) {
    const std::vector<double>& lg2 = lgam_tab(K + step + 1);
    double tk = gp_lterm(K + step, mu, lambda, lg2[K + step]);
    if (!R_finite(tk) || tk < mx - 40.0) break;
    K += step;
  }
  return std::min(K + step, cap);
}

// n x 3: log Z, mean, variance of the generalized Poisson at each mu (for
// lambda >= 0 the moments are the closed forms mu and mu/(1-lambda)^2, as in R)
// [[Rcpp::export]]
NumericMatrix gp_norm_cpp(NumericVector mu, double lambda, int cap) {
  int n = mu.size(); NumericMatrix out(n, 3);
  for (int i = 0; i < n; i++) {
    double m = mu[i];
    if (!(m > 0.0) || !R_finite(m) || lambda >= 1.0) { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; out(i, 2) = NA_REAL; continue; }
    int K = std::min(gp_kmax(m, lambda, cap), cap);
    const std::vector<double>& lg = lgam_tab(K);
    double mx = -1e300;
    std::vector<double> lt(K + 1);
    for (int k = 0; k <= K; k++) { lt[k] = gp_lterm(k, m, lambda, lg[k]); if (lt[k] > mx) mx = lt[k]; }
    double s = 0.0, s1 = 0.0, s2 = 0.0;
    for (int k = 0; k <= K; k++) { if (!R_finite(lt[k])) continue; double w = std::exp(lt[k] - mx); s += w; s1 += k * w; s2 += (double)k * k * w; }
    double e1 = s1 / s;
    out(i, 0) = mx + std::log(s);
    if (lambda >= 0.0) { out(i, 1) = m; out(i, 2) = m / ((1.0 - lambda) * (1.0 - lambda)); }
    else { out(i, 1) = e1; out(i, 2) = s2 / s - e1 * e1; }
  }
  return out;
}
