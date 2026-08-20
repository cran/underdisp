// gec_ll.cpp -- King's Generalized Event Count / Katz-family likelihood, spanning
// under-, equi-, and over-dispersion through a single dispersion parameter delta
// (the variance-to-mean ratio). Parameterize by mean mu and delta:
//   gamma = (delta - 1) / delta,   omega = mu / delta,
//   P(k+1)/P(k) = (omega + gamma k) / (k + 1).
// delta < 1 (gamma < 0): underdispersed, finite support (the CPB / binomial cell).
// delta = 1 (gamma = 0): Poisson. delta > 1 (gamma > 0): overdispersed (Katz = NB).
// The recursion gives exact moments (mean mu, variance delta*mu), which is the
// Winkelmann-Signorino-King moment correction realized directly rather than
// approximated. Computed in log-space for numerical stability.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// log P(Y=y | mu, delta); fills p0out with P(0) if non-null. Returns -1e300 when
// y lies beyond the (finite) support of an underdispersed fit.
static double gec_logpmf(int y, double mu, double delta, int max_support,
                         double* p0out, std::vector<double>& lw) {
  double gamma = (delta - 1.0) / delta;
  double omega = mu / delta;
  lw.clear(); lw.push_back(0.0);
  double mx = 0.0;
  int k = 0;
  while (k <= max_support) {
    double num = omega + gamma * (double)k;
    if (num <= 0.0) break;                          // finite support (gamma < 0)
    double next = lw[k] + std::log(num) - std::log((double)(k + 1));
    lw.push_back(next);
    if (next > mx) mx = next;
    k++;
    if (gamma >= 0.0 && next - mx < -36.0) break;   // negligible overdispersed tail
  }
  int K = (int)lw.size() - 1;
  double s = 0.0;
  for (int j = 0; j <= K; j++) s += std::exp(lw[j] - mx);
  double logZ = mx + std::log(s);
  if (p0out) *p0out = std::exp(0.0 - logZ);          // lw[0] = 0
  if (y > K) return -1e300;                          // infeasible
  return lw[y] - logZ;
}

// [[Rcpp::export]]
double gec_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset, int max_support) {
  int n = X.nrow(), p = X.ncol();
  double delta = std::exp(params[p]);
  if (delta <= 1e-8 || delta >= 1e8) return 1e10;
  double ll = 0.0;
  std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    double eta = offset[i];
    for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double lp = gec_logpmf(Y[i], std::exp(eta), delta, max_support, nullptr, lw);
    if (lp <= -1e299) return 1e10;
    ll += lp;
  }
  return -ll;
}

// [[Rcpp::export]]
NumericMatrix gec_lp0_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset, int max_support) {
  int n = X.nrow(), p = X.ncol();
  double delta = std::exp(params[p]);
  NumericMatrix out(n, 2);
  std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    double eta = offset[i];
    for (int j = 0; j < p; j++) eta += X(i, j) * params[j];
    double p0;
    out(i, 0) = gec_logpmf(Y[i], std::exp(eta), delta, max_support, &p0, lw);
    out(i, 1) = p0;
  }
  return out;
}

// n x (kmax+1) predicted-pmf matrix at (mu, delta), for scoring / calibration.
// [[Rcpp::export]]
NumericMatrix gec_pmf_cpp(NumericVector mu, double delta, int kmax, int max_support) {
  int n = mu.size();
  NumericMatrix out(n, kmax + 1);
  std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    gec_logpmf(0, mu[i], delta, max_support, nullptr, lw);
    int K = (int)lw.size() - 1;
    double mxx = -1e300; for (int j = 0; j <= K; j++) if (lw[j] > mxx) mxx = lw[j];
    double s = 0.0; for (int j = 0; j <= K; j++) s += std::exp(lw[j] - mxx);
    double logZ = mxx + std::log(s);
    for (int k = 0; k <= kmax; k++) out(i, k) = (k <= K) ? std::exp(lw[k] - logZ) : 0.0;
  }
  return out;
}

// per-observation mean of the GEC at (mu, delta), for fitted values / QoI.
// [[Rcpp::export]]
NumericVector gec_mean_cpp(NumericVector mu, double delta, int max_support) {
  int n = mu.size();
  NumericVector out(n);
  std::vector<double> lw;
  for (int i = 0; i < n; i++) {
    gec_logpmf(0, mu[i], delta, max_support, nullptr, lw);
    int K = (int)lw.size() - 1;
    double mxx = 0.0; for (int j = 0; j <= K; j++) if (lw[j] > mxx) mxx = lw[j];
    double s = 0.0, m = 0.0;
    for (int j = 0; j <= K; j++) { double w = std::exp(lw[j] - mxx); s += w; m += j * w; }
    out[i] = m / s;
  }
  return out;
}
