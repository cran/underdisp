// hp_mean.cpp -- mean-parameterized hyper-Poisson (Saez-Castillo and Conde-Sanchez 2013).
//
// The hyper-Poisson (Bardwell and Crow 1964) with rate lambda and dispersion gamma,
//   P(Y = k) = lambda^k / ((gamma)_k 1F1(1; gamma; lambda)),  (gamma)_k = Gamma(gamma + k) / Gamma(gamma),
// is underdispersed for gamma < 1, the Poisson at gamma = 1, and overdispersed for gamma > 1.
// lambda is not the mean. The mean parameterization sets mu_i = exp(x_i'beta) and recovers, for
// every observation, the rate solving E[Y | lambda_i, gamma] = mu_i. The family is an exponential
// family in log lambda with sufficient statistic k, so dE[Y]/d(log lambda) = Var(Y) > 0: the
// root is unique, and it is found by the safeguarded Newton iteration of cmp_mean.cpp.
//
// The log terms k log(lambda) - log (gamma)_k are concave in k (log (gamma)_k is convex), and
// consecutive terms satisfy t_{k+1} / t_k = lambda / (gamma + k), so they rise to the mode
// m = floor(lambda - gamma) + 1 (0 when lambda < gamma) and fall after it. The normalizer and the
// moments are summed outward from m by that ratio, each side stopping once its term falls 40 log
// units below t_m, and the moments are taken about m. A pmf whose body is wider than 1e5 terms (a
// mean of order 1e7 or more; decided before summing where a lower bound on the width already
// exceeds the cap), or whose mode does not fit an int, is refused: the root search then
// returns NA, and the likelihood treats the point as infeasible rather than summing 1e5 terms per
// Newton step. log (gamma)_k = lgamma(gamma + k) - lgamma(gamma) is evaluated as
// lgamma(k) - lbeta(gamma, k), which does not cancel when gamma is large: the overdispersed limit
// gamma -> infinity (the geometric) is where an optimizer goes on data more dispersed than it, and
// the difference of two lgamma values of order 1e10 there carried errors of 1e-6 per observation.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static inline double hp_lpoch(double gamma, double k) {
  return (k > 0.0) ? R::lgammafn(k) - R::lbeta(gamma, k) : 0.0;
}
static inline double hp_lterm(double k, double ll, double gamma) { return k * ll - hp_lpoch(gamma, k); }

static const int HP_CAP = 100000;

// log normalizer log 1F1(1; gamma; lambda), E[Y], Var(Y) at (log lambda, gamma) in one pass;
// false where the pmf is refused (see the header)
static bool hp_moments(double ll, double gamma, double& logF, double& mean, double& var) {
  const double cut = -40.0;
  double lam = std::exp(ll), md = lam - gamma;
  if (!R_finite(lam) || !(gamma > 0.0) || !(md < 1e9)) return false;
  int m = (md >= 0.0) ? (int)std::floor(md) + 1 : 0;
  if (m > 0) {
    // a body the cap cannot hold, refused before any summing: with a = gamma + m, so that
    // a - 1 <= lambda < a, the j-th log term after the mode is at least -j(j+1)/(2(a-1)) and the
    // j-th before it at least -j(j+1)/(2(a-j)), so the loops below keep at least ju and jd terms;
    // when 1 + ju + jd exceeds the cap they would refuse, and the outcome is theirs
    double a = gamma + (double)m;
    double ju = std::floor(0.5 * (std::sqrt(1.0 + 320.0 * (a - 1.0)) - 1.0)) - 1.0;   // j(j+1) <= 80(a-1)
    double jd = std::floor(0.5 * (std::sqrt(6561.0 + 320.0 * a) - 81.0)) - 1.0;       // j(j+1) <= 80(a-j)
    jd = std::min(jd, (double)m);
    if (1.0 + std::max(ju, 0.0) + std::max(jd, 0.0) > (double)HP_CAP) return false;
  }
  double s = 1.0, s1 = 0.0, s2 = 0.0, lt = 0.0;                    // sums of w, (k - m) w, (k - m)^2 w
  int body = 1;
  for (int k = m; ; k++) {                                            // upward: t_{k+1} / t_m
    lt += ll - std::log(gamma + k);
    if (lt < cut) break;
    double w = std::exp(lt), d = (double)(k + 1 - m);
    s += w; s1 += d * w; s2 += d * d * w;
    if (++body > HP_CAP) return false;
  }
  lt = 0.0;
  for (int k = m; k > 0; k--) {                                       // downward: t_{k-1} / t_m
    lt += std::log(gamma + k - 1.0) - ll;
    if (lt < cut) break;
    double w = std::exp(lt), d = (double)(k - 1 - m);
    s += w; s1 += d * w; s2 += d * d * w;
    if (++body > HP_CAP) return false;
  }
  double e1 = s1 / s;
  logF = hp_lterm((double)m, ll, gamma) + std::log(s);
  mean = (double)m + e1;
  var  = s2 / s - e1 * e1;
  return true;
}

// log lambda_i solving E[Y | lambda_i, gamma] = mu_i, one root per element of mu
// [[Rcpp::export]]
NumericVector hp_loglambda_cpp(NumericVector mu, double gamma) {
  int n = mu.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double m = mu[i];
    if (!(m > 0.0) || !R_finite(m) || !(gamma > 0.0)) { out[i] = NA_REAL; continue; }
    // start: the mean is lambda - (gamma - 1)(1 - P(0)), so lambda is near mu + gamma - 1 when
    // P(0) is small and near mu * gamma when lambda is small; for large gamma the pmf approaches
    // the geometric with ratio lambda / gamma, whose mean is mu at lambda = gamma mu / (1 + mu)
    double base = std::min(std::max(m + gamma - 1.0, m * gamma), gamma * m / (1.0 + m) + m);
    if (base < 1e-8) base = 1e-8;
    double ll = std::log(base);
    double lo = -HUGE_VAL, hi = HUGE_VAL;
    double logF, mean, var, f = HUGE_VAL;
    int refused = 0;
    for (int it = 0; it < 200; it++) {
      if (!hp_moments(ll, gamma, logF, mean, var)) {         // a refused trial rate lies above the root's
        hi = ll; f = HUGE_VAL;                               // representable range: back off toward lo,
        if (++refused > 8) break;                            // and give up on a root beyond it
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

// n x 3 matrix of (log normalizer, mean, variance) at each log lambda
// [[Rcpp::export]]
NumericMatrix hp_norm_cpp(NumericVector loglambda, double gamma) {
  int n = loglambda.size();
  NumericMatrix out(n, 3);
  double logF, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i]) || !(gamma > 0.0)) { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; out(i, 2) = NA_REAL; continue; }
    if (!hp_moments(loglambda[i], gamma, logF, mean, var)) { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; out(i, 2) = NA_REAL; continue; }
    out(i, 0) = logF; out(i, 1) = mean; out(i, 2) = var;
  }
  return out;
}
