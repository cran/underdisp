// ewp_norm.cpp -- the three-parameter exponentially weighted Poisson (Ridout and Besbeas 2004).
//
//   P(Y = k) = w_k e^{-lambda} lambda^k / k! / W(lambda, beta1, beta2),
//   w_k = exp(-beta1 (lambda - k)) for k <= lambda,  exp(-beta2 (k - lambda)) for k > lambda,
//
// a Poisson reweighted by a two-sided exponential in the distance from the centrality parameter
// lambda: beta1 and beta2 > 0 thin the lower and upper tails (underdispersion), negative values
// thicken them, and beta1 = beta2 = 0 is the Poisson. lambda is not the mean.
//
// The normalizer W is the exact sum over the support, computed on the log scale. The log terms
// are the Poisson's (concave in k) plus a piecewise-linear weight; when beta1 + beta2 < 0 the
// weight's kink at lambda is convex, so the sum starts from the larger of the Poisson's and the
// weight's reach and runs until a term falls 40 log units below the running peak, with a hard
// cap of 1e5 terms: the rule of cmp_mean.cpp and R/count_families.R. A sum the cap stops before
// its last term has fallen 40 log units below the peak is refused (NA), as there.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static inline double ewp_lterm(double k, double ll, double lam, double b1, double b2) {
  double lw = (k <= lam) ? -b1 * (lam - k) : -b2 * (k - lam);
  return k * ll - R::lgammafn(k + 1.0) - lam + lw;
}

static int ewp_kmax(double ll, double b1, double b2) {
  double lam = std::exp(ll);
  // the upper tail decays like the Poisson's times exp(-beta2 (k - lambda)): start past the
  // mode of the reweighted Poisson, lambda * exp(-beta2), and past lambda
  double reach = std::max(lam, lam * std::exp(-b2));
  if (!R_finite(reach) || reach > 1e5) reach = 1e5;
  int K = (int)std::ceil(reach + 12.0 * std::sqrt(reach + 1.0) + 30.0);
  double peak = -1e300;
  for (int k = 0; k <= std::min(K, 100000); k++) {
    double t = ewp_lterm((double)k, ll, lam, b1, b2);
    if (t > peak) peak = t;
  }
  int step = std::max(16, (int)std::min(std::ceil(reach * 0.25) + 1.0, 10000.0));
  while (K < 100000) {
    double tk = ewp_lterm((double)(K + step), ll, lam, b1, b2);
    if (!R_finite(tk) || tk < peak - 40.0) break;
    if (tk > peak) peak = tk;
    K += step;
  }
  return std::max(50, std::min(K + step, 100000));
}

// log W, E[Y], Var(Y) at (log lambda, beta1, beta2) in one pass over the support; false where the
// cap truncates the sum
static bool ewp_moments(double ll, double b1, double b2, double& logW, double& mean, double& var) {
  double lam = std::exp(ll);
  int K = ewp_kmax(ll, b1, b2);
  std::vector<double> lw(K + 1);
  double mx = -1e300;
  for (int k = 0; k <= K; k++) {
    lw[k] = ewp_lterm((double)k, ll, lam, b1, b2);
    if (lw[k] > mx) mx = lw[k];
  }
  if (K >= 100000 && (lw[K] >= mx - 40.0 || lw[K] > lw[K - 1])) return false;   // cut while large, or still rising
  double s = 0.0, s1 = 0.0, s2 = 0.0;
  for (int k = 0; k <= K; k++) {
    double w = std::exp(lw[k] - mx);
    s += w; s1 += (double)k * w;
  }
  logW = mx + std::log(s);
  mean = s1 / s;
  // the variance about the mean in a second pass: E[Y^2] - E[Y]^2 cancels when the distribution is concentrated
  for (int k = 0; k <= K; k++) { double d = (double)k - mean; s2 += d * d * std::exp(lw[k] - mx); }
  var  = s2 / s;
  return true;
}

// n x 3 matrix of (log W, mean, variance) at each log lambda
// [[Rcpp::export]]
NumericMatrix ewp_norm_cpp(NumericVector loglambda, double beta1, double beta2) {
  int n = loglambda.size();
  NumericMatrix out(n, 3);
  double logW, mean, var;
  for (int i = 0; i < n; i++) {
    if (!R_finite(loglambda[i]) || !R_finite(beta1) || !R_finite(beta2) ||
        !ewp_moments(loglambda[i], beta1, beta2, logW, mean, var)) {
      out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; out(i, 2) = NA_REAL; continue;
    }
    out(i, 0) = logW; out(i, 1) = mean; out(i, 2) = var;
  }
  return out;
}
