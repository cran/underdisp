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
// units below the peak, with a hard cap of 1e5 terms. The root search refuses a
// trial rate whose sum the cap would truncate (a mean the cap cannot represent):
// it backs off toward the bracket's lower end, gives up after eight refusals, and
// returns NA for a root it cannot solve, so the likelihood treats the point as
// infeasible. Before this rule an optimizer step to such a mean ran 200 Newton
// steps of 1e5 terms per observation (a mean-parameterized fit to zero-heavy
// overdispersed data ran for more than half an hour), at a truncated, and so
// understated, normalizer. The normalizer and the moments at a given rate
// (cmp_logZ_cpp, cmp_moments_cpp) are NA where the cap truncates the sum.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// log(k!) for k = 0..K, tabulated once and extended on demand: the sums below
// are dominated by these evaluations, and the table is a cumulative sum of logs (relative error about 3e-14 per term, 3e-9 at k = 1e5).
static const std::vector<double>& lgam_table(int K) {
  static std::vector<double> lg(1, 0.0);
  if ((int)lg.size() <= K) {
    size_t old = lg.size(); lg.resize((size_t)K + 1);
    for (size_t k = old; k <= (size_t)K; k++) lg[k] = lg[k - 1] + std::log((double)k);
  }
  return lg;
}

static int cmp_kmax(double ll, double nu, bool& capped) {
  double mode = std::exp(ll / std::max(nu, 0.05));
  capped = !R_finite(mode) || mode > 1e5;
  if (!R_finite(mode) || mode > 1e5) mode = 1e5;
  if (mode < 1.0) mode = 1.0;
  double peak = mode * ll - nu * R::lgammafn(mode + 1.0);
  if (!R_finite(peak)) peak = 0.0;
  int K = (int)std::ceil(mode);
  int step = std::max(16, (int)std::min(std::ceil(mode * 0.25) + 1.0, 10000.0));
  bool tail = false;
  while (K < 100000) {
    double tk = (double)(K + step) * ll - nu * R::lgammafn((double)(K + step) + 1.0);
    if (!R_finite(tk) || tk < peak - 40.0) { tail = true; break; }
    K += step;
  }
  if (!tail) capped = true;
  return std::max(50, std::min(K + step, 100000));
}

// log Z, E[Y], Var(Y) at (log lambda, nu) in one pass over the support; false
// where the cap truncated the sum (the root search, `stop_if_capped`, then skips
// the sum it would discard)
static bool cmp_moments(double ll, double nu, double& logZ, double& mean, double& var, bool stop_if_capped = false) {
  bool capped;
  int K = cmp_kmax(ll, nu, capped);
  if (capped && stop_if_capped) return false;
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
    s += w; s1 += (double)k * w;
  }
  logZ = mx + std::log(s);
  mean = s1 / s;
  // the variance about the mean in a second pass (a one-pass E[Y^2] - E[Y]^2 cancels at large nu)
  for (int k = 0; k <= K; k++) { double d = (double)k - mean; s2 += d * d * std::exp(lw[k] - mx); }
  var  = s2 / s;
  return !capped;
}

// log lambda_i solving E[Y | lambda_i, nu] = mu_i, one root per element of mu.
// [[Rcpp::export]]
NumericVector cmp_loglambda_cpp(NumericVector mu, double nu) {
  int n = mu.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double m = mu[i];
    if (!(m > 0.0) || !R_finite(m)) { out[i] = NA_REAL; continue; }
    // Shmueli et al. (2005) closed-form approximation as the starting value for nu >= 1. Below 1 its correction
    // (nu - 1)/(2 nu) diverges, the start falls to lambda near 1, and the first Newton steps sum tens of thousands of
    // terms (more than the cap as nu -> 0): the start interpolates instead, on the log scale, between the Poisson
    // (nu = 1, lambda = mu) and the geometric (nu -> 0, lambda = mu / (1 + mu)).
    double ll;
    if (nu >= 1.0) {
      double base = m + (nu - 1.0) / (2.0 * nu);
      if (base < 1e-8) base = 1e-8;
      ll = nu * std::log(base);
      // a mean below 1 at large nu is a Bernoulli's, whose root is log(mu / (1 - mu)); the approximation starts
      // near nu log(1/2) there, more Newton steps away than the clipped steps cover
      if (m < 1.0) ll = std::max(ll, std::log(m / (1.0 - m)));
    } else ll = nu * std::log(m) + (1.0 - nu) * std::log(m / (1.0 + m));
    double lo = -HUGE_VAL, hi = HUGE_VAL;
    double logZ, mean, var, f = HUGE_VAL;
    int refused = 0;
    for (int it = 0; it < 200; it++) {
      if (!cmp_moments(ll, nu, logZ, mean, var, true)) {     // the cap truncates this trial rate's sum:
        hi = ll; f = HUGE_VAL;                               // back off toward lo, and give up on a
        if (++refused > 8) break;                            // root beyond the representable range
        ll = R_finite(lo) ? 0.5 * (lo + hi) : ll - 1.0;
        continue;
      }
      f = mean - m;
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
    out[i] = (std::fabs(f) <= 1e-8 * std::max(1.0, m)) ? ll : NA_REAL;   // NA: no representable root
  }
  return out;
}

// log Z(lambda, nu) for each log lambda; NA where the cap truncates the sum
// (a rate the cap cannot represent: the mean is then of order 1e5 or more)
// [[Rcpp::export]]
NumericVector cmp_logZ_cpp(NumericVector loglambda, double nu) {
  int n = loglambda.size();
  NumericVector out(n);
  double logZ, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i])) { out[i] = NA_REAL; continue; }
    out[i] = cmp_moments(loglambda[i], nu, logZ, mean, var, true) ? logZ : NA_REAL;
  }
  return out;
}

// n x 2 matrix of (mean, variance) at each log lambda; NA where the cap truncates the sum
// [[Rcpp::export]]
NumericMatrix cmp_moments_cpp(NumericVector loglambda, double nu) {
  int n = loglambda.size();
  NumericMatrix out(n, 2);
  double logZ, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i]) || !cmp_moments(loglambda[i], nu, logZ, mean, var, true)) {
      out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; continue;
    }
    out(i, 0) = mean; out(i, 1) = var;
  }
  return out;
}
