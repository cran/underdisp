// cpb_pmf.h -- the CPB pmf and its derivatives, normalized by a recursion
// outward from the mode of its terms. Shared by cpb_ll.cpp (the pooled
// likelihoods) and cpb_fe.cpp (the concentrated fixed-effects likelihood and
// its gradient); the R reference implementation (.cpb_pmf_core in
// R/calibration.R) sums the lgamma form directly and is the independent check.
//
// With n = lambda/(1-alpha) and K = floor(n) the pmf is P(Y = k) = t_k / Z on
// k = 0..K, where t_k = Gamma(n+1)/(Gamma(k+1) Gamma(n-k+1)) (1-alpha)^k alpha^(n-k)
// and Z = sum_k t_k. Consecutive terms satisfy
//     t_{k+1} / t_k = (n - k)/(k + 1) * (1-alpha)/alpha,
// which is at least one exactly while k <= (n+1)(1-alpha) - 1, so the terms
// rise to the mode m = min(K, floor((n+1)(1-alpha))) and fall after it. Z is
// therefore summed outward from m by the ratio, each side stopping once its
// running term is 36 log units (a factor 2e-16) below t_m: the cost of a
// likelihood evaluation is the width of the pmf's body, not its ceiling. The
// digamma values the derivatives need follow the same recursion,
// psi(x + 1) = psi(x) + 1/x, from one evaluation at the mode.
#ifndef UNDERDISP_CPB_PMF_H
#define UNDERDISP_CPB_PMF_H
#include <Rcpp.h>
#include <cmath>

// Normalizer and (optionally) the pmf moments E[k] and E[psi(n-k+1)]:
// returns false where the point is infeasible (y above its ceiling, or the
// ceiling above max_support). lgn1 = lgamma(n+1) is returned for the callers'
// direct term.
static inline bool cpb_norm(int y, double lam, double alpha, double la, double l1a, int max_support, bool moments,
                            double& n, int& K, double& lgn1, double& logZ, double& ek, double& epsi) {
  n = lam / (1.0 - alpha);
  K = (int)std::floor(n + 1e-9);                     // integer ceilings are feasible (tolerance shared with R)
  if (y > K || K > max_support) return false;
  int m = (int)std::floor((n + 1.0) * (1.0 - alpha)); if (m > K) m = K; if (m < 0) m = 0;
  lgn1 = R::lgammafn(n + 1.0);
  double logtm = lgn1 - R::lgammafn((double)m + 1.0) - R::lgammafn(n - m + 1.0) + m * l1a + (n - m) * la;
  double r0 = (1.0 - alpha) / alpha;
  double psim = moments ? R::digamma(n - m + 1.0) : 0.0;
  double S = 1.0; ek = (double)m; epsi = psim;
  const double cut = 2.3195228323e-16;               // exp(-36)
  double t = 1.0, ps = psim;
  for (int k = m; k < K; k++) {                      // upward: t_{k+1}/t_m
    t *= (n - k) / (k + 1.0) * r0;
    if (t < cut) break;
    S += t;
    if (moments) { ps -= 1.0 / (n - k); ek += t * (k + 1.0); epsi += t * ps; }   // psi(n-k) = psi(n-k+1) - 1/(n-k)
  }
  t = 1.0; ps = psim;
  for (int k = m; k > 0; k--) {                      // downward: t_{k-1}/t_m
    t *= (double)k / ((n - k + 1.0) * r0);
    if (t < cut) break;
    S += t;
    if (moments) { ps += 1.0 / (n - k + 1.0); ek += t * (k - 1.0); epsi += t * ps; }   // psi(n-k+2) = psi(n-k+1) + 1/(n-k+1)
  }
  logZ = logtm + std::log(S);
  if (moments) { ek /= S; epsi /= S; }
  return true;
}

// log P(Y = y | lambda, alpha), zero-truncated when `truncated`; -1e300 where
// infeasible (or where the zero-truncation would condition on an empty event).
// p0 receives P(Y = 0) of the untruncated pmf.
static inline double cpb_logpmf_p0(int y, double lam, double alpha, double la, double l1a, int max_support,
                                   bool truncated, double& p0) {
  double n, lgn1, logZ, ek, epsi; int K;
  if (!cpb_norm(y, lam, alpha, la, l1a, max_support, false, n, K, lgn1, logZ, ek, epsi)) { p0 = 0.0; return -1e300; }
  double lp = lgn1 - R::lgammafn((double)y + 1.0) - R::lgammafn(n - y + 1.0) + y * l1a + (n - y) * la - logZ;
  p0 = std::exp(n * la - logZ);                      // log t_0 = n log(alpha)
  if (truncated) {
    if (p0 >= 1.0 - 1e-15) return -1e300;
    lp -= std::log1p(-p0);
  }
  return lp;
}
static inline double cpb_logpmf(int y, double lam, double alpha, double la, double l1a, int max_support, bool truncated) {
  double p0; return cpb_logpmf_p0(y, lam, alpha, la, l1a, max_support, truncated, p0);
}

// d log P(y | lambda, alpha) / d eta (eta = log lambda) and / d alpha.
// Writing psi_k = digamma(n - k + 1) and E[.] for the expectation under the
// pmf, the terms common to every k cancel and
//   d/deta   = n (E[psi_k] - psi_y),
//   d/dalpha = -(y - E[k]) (1/(1-alpha) + 1/alpha) - n/(1-alpha) (psi_y - E[psi_k]);
// with zero truncation the same quantities at y = 0 enter through
//   -d log(1 - p0) = p0/(1-p0) d log p0.
static inline bool cpb_dlogpmf(int y, double lam, double alpha, double la, double l1a, int max_support, bool truncated,
                               double& deta, double& dalpha) {
  double n, lgn1, logZ, ek, epsi; int K;
  if (!cpb_norm(y, lam, alpha, la, l1a, max_support, true, n, K, lgn1, logZ, ek, epsi)) return false;
  double psy = R::digamma(n - y + 1.0);
  double c = 1.0 / (1.0 - alpha) + 1.0 / alpha, r = n / (1.0 - alpha);
  deta   = n * (epsi - psy);
  dalpha = -((double)y - ek) * c - r * (psy - epsi);
  if (truncated) {
    double p0 = std::exp(n * la - logZ);
    if (p0 >= 1.0 - 1e-15) return false;
    double f = p0 / (1.0 - p0), ps0 = R::digamma(n + 1.0);
    deta   += f * (n * (epsi - ps0));
    dalpha += f * (-(0.0 - ek) * c - r * (ps0 - epsi));
  }
  return true;
}
#endif
