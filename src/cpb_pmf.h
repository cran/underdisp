// cpb_pmf.h -- the CPB pmf and its derivatives, normalized by a recursion
// outward from the mode of its terms. Shared by cpb_ll.cpp (the pooled
// likelihoods) and cpb_fe.cpp (the concentrated fixed-effects likelihood and
// its gradient); the R reference implementation (.cpb_pmf_core in
// R/calibration.R) is the independent check.
//
// Two supports. King's CPB (King 1989; Winkelmann, Signorino and King 1995),
// the default, runs over k = 0..ceiling(n), the points at which the Katz
// recursion's numerator lambda - k(1-alpha) is still positive; it is the GEC
// of gec_ll.cpp at delta = alpha < 1, and its likelihood is continuous in the
// parameters (the entering top point's mass grows from zero). The "floor"
// support, k = 0..floor(n), is the package's earlier definition (0.1.x), kept
// unchanged; its likelihood jumps where a ceiling crosses an integer.
//
// Floor support. With n = lambda/(1-alpha) and K = floor(n) the pmf is P(Y = k) = t_k / Z on
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
#include "support_cap.h"

// Normalizer and (optionally) the pmf moments E[k] and E[psi(n-k+1)]:
// returns false where the point is infeasible (y above its ceiling, or the
// ceiling above max_support). lgn1 = lgamma(n+1) is returned for the callers'
// direct term.
static inline bool cpb_norm(int y, double lam, double alpha, double la, double l1a, int max_support, bool moments,
                            double& n, int& K, double& lgn1, double& logZ, double& ek, double& epsi) {
  n = lam / (1.0 - alpha);
  K = cap_support_int(std::floor(n + 1e-9), max_support);   // integer ceilings are feasible (tolerance shared with R)
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

// ---- King's support, 0..ceiling(n) -------------------------------------------
// The terms are carried relative to the mode m and stepped by
//     t_{k+1} / t_k = (lambda - k(1-alpha)) / (alpha (k+1)),
// whose numerator is the Katz recursion's; the support is every k at which the
// numerator before it is positive, i.e. k = 0..ceiling(n - 1e-9) (the
// tolerance is the floor support's, so both treat an integer n alike). Nothing
// here forms n or Gamma(n + 1), so the pmf stays exact as alpha approaches one
// and n diverges: the Poisson limit is approached continuously rather than
// refused, and max_support only bounds the width of the pmf's body. Each side
// is summed until its term falls 36 log units below t_m. The terms fall
// monotonically away from the mode, so inside that body they lie in
// [exp(-36), 1] and are stepped as a running product of the ratios, with no log
// or exp per term (as on the floor support); past it the walk is continued in
// log space when the observed count lies beyond it, for at most TAIL_WALK
// further steps; a count farther out (or a zero the walk does not reach) takes
// its term in closed form,
//     log(t_k / t_m) = lgamma(m+1) + lgamma(n-m+1) - lgamma(k+1) - lgamma(n-k+1)
//                      + (k - m) log((1-alpha)/alpha),
// with each lgamma difference evaluated as a log Pochhammer symbol,
// lgamma(a + d) - lgamma(a) = lgamma(d) - lbeta(a, d), which does not cancel when
// n is large (as alpha approaches one, n reaches 1e12 and more, and the plain
// differences of lgamma values of that size were wrong by whole log units; a
// zero-inflated fit then saw P(0) > 1 and a positive log-likelihood). alpha
// within 1e-12 of one is infeasible here as in every caller, and a normalized
// log probability above zero is refused as a numerical failure.
// Without the bound a trial point with a large rate and a tiny alpha (a body a
// few points wide near n) walked from the mode to a small count term by term:
// 1e8 steps per observation. The digamma differences the derivatives need,
// psi_k = digamma(n - k + 1), are carried relative to psi_m by
// psi(x + 1) = psi(x) + 1/x, with 1/(n - k) = (1-alpha)/(lambda - k(1-alpha)),
// and in closed form beside the closed-form terms.
struct CpbCeil {
  double ly;      // log P(Y = y)
  double l0;      // log P(Y = 0)
  double n;       // lambda / (1 - alpha)
  double ek;      // E[k]
  double dpy;     // E[psi_k] - psi_y
  double dp0;     // E[psi_k] - psi_0
};
static const int TAIL_WALK = 500;
// log Gamma(a + d) / Gamma(a), d >= 0
static inline double cpb_lpoch(double a, double d) { return (d > 0.0) ? R::lgammafn(d) - R::lbeta(a, d) : 0.0; }
// digamma(a + d) - digamma(a), d >= 0: the asymptotic series once a is large
static inline double cpb_psidiff(double a, double d) {
  if (!(d > 0.0)) return 0.0;
  if (a < 1e6) return R::digamma(a + d) - R::digamma(a);
  double b = a + d;
  return std::log1p(d / a) + 0.5 / a - 0.5 / b + 1.0 / (12.0 * a * a) - 1.0 / (12.0 * b * b);
}
static inline bool cpb_norm_ceil(int y, double lam, double alpha, int max_support, bool moments, CpbCeil& r) {
  if (!(lam > 0.0) || !R_finite(lam) || !(alpha > 1e-12) || !(alpha < 1.0 - 1e-12)) return false;
  const double om = 1.0 - alpha, tol = 1e-9 * om, ecut = 2.3195228302435696e-16;   // exp(-36)
  if (y < 0 || (y > 0 && !(lam - (y - 1.0) * om > tol))) return false;   // y above ceiling(n)
  if (!(lam < 1e9)) return false;                                     // the mode must fit an int
  int m = (lam > alpha) ? (int)std::floor(lam - alpha) + 1 : 0;        // the terms rise while k <= lambda - alpha
  while (m > 0 && !(lam - (m - 1.0) * om > tol)) m--;                  // keep the mode inside the support
  double S = 1.0, ek = (double)m, eps = 0.0, ly = 0.0, l0 = 0.0, psy = 0.0, ps0 = 0.0;
  bool got_y = (y == m), got_0 = (m == 0);
  int body = 1, out = 0;
  double t = 1.0, lt = 0.0, ps = 0.0;                                  // t in the body, lt = log t past it
  bool tail = false;
  for (int k = m; ; k++) {                                             // upward
    double num = lam - k * om;
    if (!(num > tol)) break;                                           // k + 1 lies outside the support
    double q = num / (alpha * (k + 1.0));
    if (moments) ps -= om / num;                                       // psi_{k+1} = psi_k - 1/(n - k)
    if (!tail) {
      double t1 = t * q;
      if (t1 >= ecut) {
        t = t1; S += t;
        if (k + 1 == y) { ly = std::log(t); psy = ps; got_y = true; }
        if (moments) { ek += t * (k + 1.0); eps += t * ps; }
        if (++body > max_support) return false;                        // a body wider than the guard
        continue;
      }
      lt = std::log(t) + std::log(q); tail = true;
    } else lt += std::log(q);
    if (k + 1 == y) { ly = lt; psy = ps; got_y = true; }
    if (k + 1 >= y || ++out > TAIL_WALK) break;
  }
  t = 1.0; lt = 0.0; ps = 0.0; out = 0; tail = false;
  for (int k = m; k > 0; k--) {                                        // downward
    double den = lam - (k - 1.0) * om;                                 // positive: k - 1 < m lies in the support
    double q = alpha * k / den;
    if (moments) ps += om / den;                                       // psi_{k-1} = psi_k + 1/(n - k + 1)
    if (!tail) {
      double t1 = t * q;
      if (t1 >= ecut) {
        t = t1; S += t;
        if (k - 1 == y) { ly = std::log(t); psy = ps; got_y = true; }
        if (k - 1 == 0) { l0 = std::log(t); ps0 = ps; got_0 = true; }
        if (moments) { ek += t * (k - 1.0); eps += t * ps; }
        if (++body > max_support) return false;
        continue;
      }
      lt = std::log(t) + std::log(q); tail = true;
    } else lt += std::log(q);
    if (k - 1 == y) { ly = lt; psy = ps; got_y = true; }
    if (k - 1 == 0) { l0 = lt; ps0 = ps; got_0 = true; }
    if (k - 1 <= y || ++out > TAIL_WALK) break;
  }
  r.n = lam / om;
  if (!got_y || !got_0) {                                              // terms the walks left: closed form
    const double lr = std::log(om / alpha), nn = r.n;
    auto lrel = [&](int k) {                                           // log(t_k / t_m)
      double d = std::fabs((double)k - (double)m);
      return (k < m) ? cpb_lpoch(k + 1.0, d) - cpb_lpoch(nn - m + 1.0, d) - d * lr
                     : -cpb_lpoch(m + 1.0, d) + cpb_lpoch(nn - k + 1.0, d) + d * lr;
    };
    auto prel = [&](int k) {                                           // psi_k - psi_m
      double d = std::fabs((double)k - (double)m);
      return (k < m) ? cpb_psidiff(nn - m + 1.0, d) : -cpb_psidiff(nn - k + 1.0, d);
    };
    if (!got_y) { ly = lrel(y); if (moments) psy = prel(y); }
    if (!got_0) { l0 = lrel(0); if (moments) ps0 = prel(0); }          // P(0) below exp(-36) of the mode
  }
  double lS = std::log(S);
  r.ly = ly - lS; r.l0 = l0 - lS;
  if (!(r.ly <= 1e-10) || !(r.l0 <= 1e-10)) return false;             // a log probability above zero
  r.ek = ek / S;
  double epsr = eps / S;                                               // E[psi_k] - psi_m
  r.dpy = epsr - psy; r.dp0 = epsr - ps0;
  return true;
}

// log P(Y = y), zero-truncated when `truncated`; -1e300 where infeasible. p0
// receives P(Y = 0) of the untruncated pmf.
static inline double cpb_logpmf_p0_ceil(int y, double lam, double alpha, int max_support, bool truncated, double& p0) {
  CpbCeil r;
  if (!cpb_norm_ceil(y, lam, alpha, max_support, false, r)) { p0 = 0.0; return -1e300; }
  p0 = std::exp(r.l0);
  double lp = r.ly;
  if (truncated) {
    if (p0 >= 1.0 - 1e-15) return -1e300;
    lp -= std::log1p(-p0);
  }
  return lp;
}

// d log P(y) / d eta and / d alpha, the formulas of cpb_dlogpmf on King's
// support (they hold between the integer crossings, where the support is fixed;
// the likelihood is continuous across a crossing and its slope jumps there).
static inline bool cpb_dlogpmf_ceil(int y, double lam, double alpha, int max_support, bool truncated,
                                    double& deta, double& dalpha) {
  CpbCeil r;
  if (!cpb_norm_ceil(y, lam, alpha, max_support, true, r)) return false;
  double om = 1.0 - alpha, c = 1.0 / om + 1.0 / alpha, rr = r.n / om;
  deta   = r.n * r.dpy;
  dalpha = -((double)y - r.ek) * c + rr * r.dpy;
  if (truncated) {
    double p0 = std::exp(r.l0);
    if (p0 >= 1.0 - 1e-15) return false;
    double f = p0 / (1.0 - p0);
    deta   += f * (r.n * r.dp0);
    dalpha += f * (r.ek * c + rr * r.dp0);
  }
  return true;
}

// ---- dispatch on the support (floor_support = the 0.1.x definition) -----------
static inline double cpb_logpmf_p0_s(bool floor_support, int y, double lam, double alpha, double la, double l1a,
                                     int max_support, bool truncated, double& p0) {
  return floor_support ? cpb_logpmf_p0(y, lam, alpha, la, l1a, max_support, truncated, p0)
                       : cpb_logpmf_p0_ceil(y, lam, alpha, max_support, truncated, p0);
}
static inline double cpb_logpmf_s(bool floor_support, int y, double lam, double alpha, double la, double l1a,
                                  int max_support, bool truncated) {
  double p0; return cpb_logpmf_p0_s(floor_support, y, lam, alpha, la, l1a, max_support, truncated, p0);
}
static inline bool cpb_dlogpmf_s(bool floor_support, int y, double lam, double alpha, double la, double l1a,
                                 int max_support, bool truncated, double& deta, double& dalpha) {
  return floor_support ? cpb_dlogpmf(y, lam, alpha, la, l1a, max_support, truncated, deta, dalpha)
                       : cpb_dlogpmf_ceil(y, lam, alpha, max_support, truncated, deta, dalpha);
}
#endif
