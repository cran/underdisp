// cmp_mean.cpp -- mean-parameterized Conway-Maxwell-Poisson (Huang 2017).
//
// The standard COM-Poisson is indexed by a rate lambda and a dispersion nu,
//   P(Y = k) = lambda^k / (k!)^nu / Z(lambda, nu),  Z = sum_j lambda^j / (j!)^nu,
// and lambda is not the mean. Huang's mean parameterization sets the mean mu_i
// = exp(x_i'beta) and recovers, for every observation, the rate lambda_i that
// solves E[Y | lambda_i, nu] = mu_i. Because dE[Y]/d(log lambda) = Var(Y) > 0,
// the mean is strictly increasing in log lambda and the root is unique; it is
// found by a safeguarded Newton iteration on log lambda.
//
// The normalizing sum is truncated by the same convergence rule as the R-side
// rate parameterization (R/compois.R): the log terms k*log(lambda) - nu*lgamma(k+1)
// are unimodal in k, so the sum runs past the mode until a term falls 40 log
// units below the peak, with a hard cap of 1e5 terms.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// log(k!) for k = 0..K, tabulated once and extended on demand: the sums below
// are dominated by these evaluations, and the table is exact (cumulative log).
static const std::vector<double>& lgam_table(int K) {
  static std::vector<double> lg(1, 0.0);
  if ((int)lg.size() <= K) {
    size_t old = lg.size(); lg.resize((size_t)K + 1);
    for (size_t k = old; k <= (size_t)K; k++) lg[k] = lg[k - 1] + std::log((double)k);
  }
  return lg;
}

static int cmp_kmax(double ll, double nu) {
  double mode = std::exp(ll / std::max(nu, 0.05));
  if (!R_finite(mode) || mode > 1e5) mode = 1e5;
  if (mode < 1.0) mode = 1.0;
  double peak = mode * ll - nu * R::lgammafn(mode + 1.0);
  if (!R_finite(peak)) peak = 0.0;
  int K = (int)std::ceil(mode);
  int step = std::max(16, (int)std::min(std::ceil(mode * 0.25) + 1.0, 10000.0));
  while (K < 100000) {
    double tk = (double)(K + step) * ll - nu * R::lgammafn((double)(K + step) + 1.0);
    if (!R_finite(tk) || tk < peak - 40.0) break;
    K += step;
  }
  return std::max(50, std::min(K + step, 100000));
}

// log Z, E[Y], Var(Y) at (log lambda, nu) in one pass over the support
static void cmp_moments(double ll, double nu, double& logZ, double& mean, double& var) {
  int K = cmp_kmax(ll, nu);
  const std::vector<double>& lg = lgam_table(K);
  std::vector<double> lw(K + 1);
  double mx = -1e300;
  for (int k = 0; k <= K; k++) {
    lw[k] = (double)k * ll - nu * lg[k];
    if (lw[k] > mx) mx = lw[k];
  }
  double s = 0.0, s1 = 0.0, s2 = 0.0;
  for (int k = 0; k <= K; k++) {
    double w = std::exp(lw[k] - mx);
    s += w; s1 += (double)k * w; s2 += (double)k * (double)k * w;
  }
  logZ = mx + std::log(s);
  mean = s1 / s;
  var  = s2 / s - mean * mean;
}

// log lambda_i solving E[Y | lambda_i, nu] = mu_i, one root per element of mu.
// [[Rcpp::export]]
NumericVector cmp_loglambda_cpp(NumericVector mu, double nu) {
  int n = mu.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double m = mu[i];
    if (!(m > 0.0) || !R_finite(m)) { out[i] = NA_REAL; continue; }
    // Shmueli et al. (2005) closed-form approximation as the starting value
    double base = m + (nu - 1.0) / (2.0 * nu);
    if (base < 1e-8) base = 1e-8;
    double ll = nu * std::log(base);
    double lo = -HUGE_VAL, hi = HUGE_VAL;
    double logZ, mean, var;
    for (int it = 0; it < 200; it++) {
      cmp_moments(ll, nu, logZ, mean, var);
      double f = mean - m;
      if (std::fabs(f) <= 1e-11 * std::max(1.0, m)) break;
      if (f > 0.0) hi = ll; else lo = ll;                    // mean is increasing in ll
      double step = f / std::max(var, 1e-12);
      if (std::fabs(step) > 5.0) step = (step > 0.0) ? 5.0 : -5.0;
      double nl = ll - step;
      if (!(nl > lo && nl < hi) || !R_finite(nl)) {          // safeguard: bisect when bracketed
        if (R_finite(lo) && R_finite(hi)) nl = 0.5 * (lo + hi);
        else nl = ll - ((f > 0.0) ? 1.0 : -1.0);
      }
      ll = nl;
    }
    out[i] = ll;
  }
  return out;
}

// log Z(lambda, nu) for each log lambda
// [[Rcpp::export]]
NumericVector cmp_logZ_cpp(NumericVector loglambda, double nu) {
  int n = loglambda.size();
  NumericVector out(n);
  double logZ, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i])) { out[i] = NA_REAL; continue; }
    cmp_moments(loglambda[i], nu, logZ, mean, var);
    out[i] = logZ;
  }
  return out;
}

// n x 2 matrix of (mean, variance) at each log lambda
// [[Rcpp::export]]
NumericMatrix cmp_moments_cpp(NumericVector loglambda, double nu) {
  int n = loglambda.size();
  NumericMatrix out(n, 2);
  double logZ, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i])) { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; continue; }
    cmp_moments(loglambda[i], nu, logZ, mean, var);
    out(i, 0) = mean; out(i, 1) = var;
  }
  return out;
}
