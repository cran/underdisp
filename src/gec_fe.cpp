// gec_fe.cpp -- King GEC / Katz-family regression with high-dimensional unit
// fixed effects, concentrated (profiled) out of the likelihood exactly as in
// cpb_fe.cpp; the outer optimizer sees only (beta, log delta). Data SORTED BY UNIT.
//
// The inner problem. For delta < 1 the Katz support is finite,
// 0..ceil(mu_t / (1 - delta)), but unlike the CPB the likelihood does not jump
// when a ceiling crosses an integer: the entering support point's mass is
// proportional to omega + gamma k and grows from zero, so the unit objective
// is continuous in the intercept with a kink (a downward jump of the slope) at
// each crossing, and it tends to -infinity at the feasibility floor, where an
// observation's ceiling reaches y - 1. Its smooth part is unimodal and the
// kinks only lower the slope, so the objective is unimodal: golden section on
// a bracket around the unit's Poisson intercept (a warm start may be supplied;
// the package evaluates from cold), expanded while the maximum sits at an
// edge, finds the maximum, which may lie at a kink. For delta >= 1 the support
// is unbounded and the objective is smooth. The gradient of the concentrated
// objective follows the envelope theorem; at a kink the intercept tracks the
// crossing a* = log(k (1-delta)) - o_{t*}, so the total derivative adds
// (dl/da) (da*/dtheta) with da*/dbeta = -x_{t*} and
// da*/dlog(delta) = -delta/(1-delta); the combination is the same from both
// sides of the kink.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

static double gec_logpmf_fe(int y, double mu, double delta, int max_support) {
  if (!R_finite(mu) || mu <= 0.0) return -1e300;      // an overflowed trial point is infeasible
  double gamma = (delta - 1.0) / delta, omega = mu / delta;
  static std::vector<double> lw, lk;                   // work buffer and log(k+1) table, reused across calls
  if ((int)lw.size() < max_support + 2) {
    lw.resize(max_support + 2); lk.resize(max_support + 2);
    for (int j = 0; j < (int)lk.size(); j++) lk[j] = std::log((double)(j + 1));
  }
  lw[0] = 0.0;
  double mx = 0.0; int k = 0;
  while (k <= max_support) {
    double num = omega + gamma * (double)k;
    if (num <= 0.0) break;
    double next = lw[k] + std::log(num) - lk[k];
    lw[k + 1] = next; if (next > mx) mx = next; k++;
    // the terms are unimodal in k (their ratio (omega + gamma k)/(k+1) is
    // decreasing), so the sum stops 36 log units below the peak once past y
    if (k >= y && next - mx < -36.0) break;
  }
  int K = k;
  if (y > K) return -1e300;
  double s = 0.0; for (int j = 0; j <= K; j++) s += std::exp(lw[j] - mx);
  return lw[y] - (mx + std::log(s));
}

// d log P(y | mu, delta) / d eta (eta = log mu) and / d log(delta). With
// lw[k] = sum_{j<k} [log(omega + gamma j) - log(j+1)], omega = mu/delta and
// gamma = 1 - 1/delta, the partial sums c_eta[k] = sum_{j<k} omega/(omega + gamma j)
// and c_delta[k] = sum_{j<k} (-omega + j/delta)/(omega + gamma j) give
//   d/deta = c_eta[y] - E[c_eta],   d/dlog(delta) = c_delta[y] - E[c_delta].
static bool gec_dlogpmf_fe(int y, double mu, double delta, int max_support, double& deta, double& ddelta) {
  if (!R_finite(mu) || mu <= 0.0) return false;
  double gamma = (delta - 1.0) / delta, omega = mu / delta;
  static std::vector<double> lw, ce, cd, lk;
  if ((int)lw.size() < max_support + 2) {
    lw.resize(max_support + 2); ce.resize(max_support + 2); cd.resize(max_support + 2); lk.resize(max_support + 2);
    for (int j = 0; j < (int)lk.size(); j++) lk[j] = std::log((double)(j + 1));
  }
  lw[0] = 0.0; ce[0] = 0.0; cd[0] = 0.0;
  double mx = 0.0; int k = 0;
  while (k <= max_support) {
    double num = omega + gamma * (double)k;
    if (num <= 0.0) break;
    lw[k + 1] = lw[k] + std::log(num) - lk[k];
    ce[k + 1] = ce[k] + omega / num;
    cd[k + 1] = cd[k] + (-omega + (double)k / delta) / num;
    if (lw[k + 1] > mx) mx = lw[k + 1]; k++;
    if (k >= y && lw[k] - mx < -36.0) break;
  }
  int K = k;
  if (y > K) return false;
  double s = 0.0, ee = 0.0, ed = 0.0;
  for (int j = 0; j <= K; j++) { double w = std::exp(lw[j] - mx); s += w; ee += w * ce[j]; ed += w * cd[j]; }
  deta = ce[y] - ee / s; ddelta = cd[y] - ed / s;
  return true;
}

namespace {   // internal linkage: cpb_fe.cpp defines a different Unit
struct Unit {
  const IntegerVector& Y; const std::vector<double>& o; const NumericVector& w;
  int lo, hi; double delta; int max_support;
  double ll(double a) const {
    double v = 0.0;
    for (int t = lo; t < hi; t++) {
      double lp = gec_logpmf_fe(Y[t], std::exp(a + o[t]), delta, max_support);
      if (lp <= -1e299 || !R_finite(lp)) return -1e300;
      v += w[t] * lp;
    }
    return v;
  }
  // integer crossings bounding the piece that contains a (delta < 1):
  // prev <= a < next, and the observations that own them
  void tooth(double a, double& prev, double& next, int& prev_t, int& next_t) const {
    prev = -1e300; next = 1e300; prev_t = -1; next_t = -1;
    for (int t = lo; t < hi; t++) {
      double N = std::exp(a + o[t]) / (1.0 - delta);
      int K = (int)std::floor(N + 1e-9);
      if (K >= 1) { double b = std::log((double)K * (1.0 - delta)) - o[t]; if (b > prev) { prev = b; prev_t = t; } }
      double b2 = std::log((double)(K + 1) * (1.0 - delta)) - o[t]; if (b2 < next) { next = b2; next_t = t; }
    }
  }
};
}  // namespace

static const double KINK_EPS = 1e-8;
static const double GR = 0.6180339887498949;

// golden section for a maximum on [L, H]: at most `iters` iterations, stopping
// once the bracket is narrower than 1e-7 on the log scale
static double golden_max(const Unit& U, double L, double H, int iters, double& a_at) {
  double c = H - GR * (H - L), d = L + GR * (H - L), fc = U.ll(c), fd = U.ll(d);
  for (int it = 0; it < iters && (H - L) > 1e-7; it++) {
    if (fc < fd) { L = c; c = d; fc = fd; d = L + GR * (H - L); fd = U.ll(d); }
    else         { H = d; d = c; fd = fc; c = H - GR * (H - L); fc = U.ll(c); }
  }
  a_at = 0.5 * (L + H);
  return U.ll(a_at);
}

// returns the unit's maximized log-likelihood; a_star = the intercept; edge_t =
// the row whose ceiling crossing the intercept sits at (-1 for an interior optimum)
static double maximize_unit(const Unit& U, double a0, double& a_star, int& edge_t, int inner_it) {
  double alo = -1e300; int ymax = 0; double sy = 0.0, se = 0.0;
  bool finite_support = U.delta < 1.0;
  edge_t = -1;
  for (int t = U.lo; t < U.hi; t++) {
    sy += U.w[t] * (double)U.Y[t]; se += U.w[t] * std::exp(U.o[t]);
    if (U.Y[t] > 0) {
      if (U.Y[t] > ymax) ymax = U.Y[t];
      if (finite_support && U.Y[t] >= 2) {         // support 0..ceil(N): need N > y - 1
        double b = std::log((1.0 - U.delta) * (double)(U.Y[t] - 1)) - U.o[t];
        if (b > alo) alo = b;
      }
    }
  }
  if (ymax == 0) { a_star = -30.0; return 0.0; }
  double apois = (sy > 0.0 && se > 0.0) ? std::log(sy / se) : 0.0;
  double L0 = finite_support ? alo + 1e-9 : -30.0;
  double centre;
  if (R_finite(a0)) centre = std::max(a0, L0);
  else {
    double L = std::max(L0, apois - 2.0), H = std::max(L + 0.5, apois + 2.0);
    double bv = -1e300; centre = L;
    for (int ext = 0; ext < 3; ext++) {
      int G = (int)std::ceil((H - L) / 0.5);
      for (int g = 0; g <= G; g++) {
        double a = L + (H - L) * (double)g / (double)G;
        double v = U.ll(a);
        if (v > bv) { bv = v; centre = a; }
      }
      if (centre >= H - 1e-9) { L = H; H = H + 2.0; continue; }
      if (centre <= L + 1e-9 && L > L0 + 1e-9) { H = L; L = std::max(L0, L - 2.0); continue; }
      break;
    }
  }
  // unimodal objective: golden section on a bracket that expands while the
  // maximum sits at an edge
  double best_a = centre, best_v = U.ll(centre);
  double L = std::max(L0, centre - 0.5), H = centre + 0.5, ag, vg;
  for (int ext = 0; ext < 6; ext++) {
    vg = golden_max(U, L, H, inner_it, ag);
    if (vg > best_v) { best_v = vg; best_a = ag; }
    if (ag >= H - 1e-4) { L = H - 0.1; H = H + 1.0; continue; }
    if (ag <= L + 1e-4 && L > L0 + 1e-9) { H = L + 0.1; L = std::max(L0, L - 1.0); continue; }
    break;
  }
  if (best_v <= -1e299) { a_star = L0; return U.ll(L0); }
  // a maximum at a ceiling crossing (a kink) tracks that crossing as the
  // parameters move: flag it for the gradient (finite support only)
  if (finite_support) {
    double prev, next; int pt, nt; U.tooth(best_a, prev, next, pt, nt);
    if (best_a >= next - 1e-5 && nt >= 0) {
      double aH = next - KINK_EPS, vH = U.ll(aH);
      if (vH >= best_v - 1e-12) { best_v = vH; best_a = aH; }
      edge_t = nt;
    } else if (best_a <= prev + 1e-5 && pt >= 0) {
      double aL = std::max(L0, prev + KINK_EPS), vL = U.ll(aL);
      if (vL >= best_v - 1e-12) { best_v = vL; best_a = aL; }
      edge_t = pt;
    }
  }
  a_star = best_a;
  return best_v;
}

static void compute_offsets(NumericVector params, NumericMatrix X, NumericVector uoff, std::vector<double>& o) {
  int n = X.nrow(), p = X.ncol();
  for (int i = 0; i < n; i++) { double e = uoff[i]; for (int j = 0; j < p; j++) e += X(i, j) * params[j]; o[i] = e; }
}

// Concentrated negative log-likelihood at (beta, log delta); `awarm` as in
// cpb_fe_nll_cpp (length n_units = warm starts, length 0 = cold). Returns
// list(nll, a, edge).
// [[Rcpp::export]]
List gec_fe_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                    NumericVector w, IntegerVector ustart, int n_units, int max_support, int inner_it,
                    NumericVector awarm) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  bool warm = awarm.size() == n_units;
  IntegerVector edge(n_units, -1);
  double delta = std::exp(params[p]);
  if (delta <= 1e-8 || delta >= 1e8) return List::create(1e10, awarm, edge);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector a(n_units);
  double total = 0.0, astar; int et;
  for (int u = 0; u < n_units; u++) {
    if ((u & 15) == 0) Rcpp::checkUserInterrupt();
    Unit U = { Y, o, w, ustart[u], ustart[u + 1], delta, max_support };
    double a0 = warm ? awarm[u] : NA_REAL;
    double ll = maximize_unit(U, a0, astar, et, inner_it);
    if (ll <= -1e299) return List::create(1e10, awarm, edge);
    a[u] = astar; edge[u] = et; total += ll;
  }
  return List::create(-total, a, edge);
}

// Gradient of the concentrated negative log-likelihood with respect to
// (beta, log delta) at the solved intercepts `a`.
// [[Rcpp::export]]
NumericVector gec_fe_grad_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                              NumericVector w, IntegerVector ustart, int n_units, int max_support,
                              NumericVector a, IntegerVector edge) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || a.size() != n_units || edge.size() != n_units)
    Rcpp::stop("internal error: gradient inputs do not match the design");
  double delta = std::exp(params[p]);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector G(p + 1);
  for (int u = 0; u < n_units; u++) {
    double Ga = 0.0, Gd = 0.0; std::vector<double> Gb(p, 0.0);
    for (int t = ustart[u]; t < ustart[u + 1]; t++) {
      double de, dd;
      if (!gec_dlogpmf_fe(Y[t], std::exp(a[u] + o[t]), delta, max_support, de, dd))
        return NumericVector(p + 1, NA_REAL);
      Ga += w[t] * de; Gd += w[t] * dd;
      for (int j = 0; j < p; j++) Gb[j] += w[t] * de * X(t, j);
    }
    if (edge[u] >= 0 && delta < 1.0) {                    // the intercept tracks a ceiling crossing
      int ts = edge[u];
      for (int j = 0; j < p; j++) Gb[j] -= Ga * X(ts, j);
      Gd -= Ga * delta / (1.0 - delta);
    }
    for (int j = 0; j < p; j++) G[j] -= Gb[j];
    G[p] -= Gd;
  }
  return G;
}

// Unit intercepts at (beta, log delta) from a cold search for every unit.
// [[Rcpp::export]]
NumericVector gec_fe_intercepts_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                                    NumericVector w, IntegerVector ustart, int n_units, int max_support, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  double delta = std::exp(params[p]);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector fe(n_units); double astar; int et;
  for (int u = 0; u < n_units; u++) {
    Unit U = { Y, o, w, ustart[u], ustart[u + 1], delta, max_support };
    maximize_unit(U, NA_REAL, astar, et, inner_it);
    fe[u] = astar;
  }
  return fe;
}
