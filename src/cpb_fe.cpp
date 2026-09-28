// cpb_fe.cpp -- CPB regression with high-dimensional unit fixed effects via a
// CONCENTRATED (profiled) likelihood. Each unit's intercept is solved by a 1-D inner
// maximization, so the outer optimizer sees only (beta, alpha). No unit dummies.
// Data must be SORTED BY UNIT; ustart[u]..ustart[u+1]-1 are unit u's rows.
//
// The inner problem. Observation t has rate lambda_t = exp(a + o_t) and ceiling
// N_t = lambda_t / (1 - alpha), and its pmf is renormalized over 0..floor(N_t).
// As a increases through a breakpoint a = log(k (1-alpha)) - o_t (N_t = k) the
// support gains the point k, the normalizer grows, and the log pmf of the
// observed y_t < k drops by about (1-alpha)^k. The unit log-likelihood is
// therefore a saw-tooth in a: smooth between breakpoints, with a downward jump
// at every breakpoint, and its smooth part is unimodal. Its supremum over a is
// attained at the interior critical point of the tooth that contains the peak
// of the smooth part, as the left limit at a breakpoint, or at the feasibility
// floor (the left end of the first tooth, where an observation's count equals
// its ceiling). The search below (i) locates the peak of the smooth part
// (a coarse grid around the unit's Poisson intercept; a warm start may be
// supplied, but the package evaluates every unit from cold, since a warm start
// carried from a distant trial point can leave a unit on the wrong tooth),
// (ii) evaluates the
// left limit at every breakpoint within a window of that peak, and (iii)
// refines by golden section inside the tooth holding the peak and inside the
// tooth ending at the best breakpoint. Solving the inner problem to its
// supremum makes the concentrated likelihood continuous in alpha (a change in
// alpha shifts every breakpoint of a unit rigidly, so the maximizing tooth
// slides; kinks where it switches). It is NOT continuous in beta: at a unit
// whose supremum sits on its feasibility floor the intercept cannot retreat,
// so when another observation's breakpoint crosses that floor the supremum
// drops by about (1-alpha)^k. The outer optimizer treats it accordingly (a
// multistart with restarts; see .fe_outer in R/utils.R).
//
// The gradient. By the envelope theorem the gradient of the concentrated
// log-likelihood is the partial gradient at the solved intercepts; where a
// unit's supremum is the left limit at a breakpoint of observation t*, the
// intercept tracks that breakpoint, a* = log(k (1-alpha)) - o_{t*}, so the
// total derivative adds (dl/da) (da*/dtheta) with da*/dbeta = -x_{t*} and
// da*/dalpha = -1/(1-alpha).
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "cpb_pmf.h"
using namespace Rcpp;

// Everything the inner search needs about one unit. Internal linkage: gec_fe.cpp
// defines a different Unit, and two definitions of one class name with external
// linkage would share their inline member functions across the files.
namespace {
struct Unit {
  const IntegerVector& Y; const std::vector<double>& o; const NumericVector& w;
  int lo, hi; double alpha, la, l1a; int max_support; bool truncated;
  double ll(double a) const {
    double v = 0.0;
    for (int t = lo; t < hi; t++) {
      double lp = cpb_logpmf(Y[t], std::exp(a + o[t]), alpha, la, l1a, max_support, truncated);
      if (lp <= -1e299 || !R_finite(lp)) return -1e300;
      v += w[t] * lp;
    }
    return v;
  }
  // breakpoints bounding the tooth that contains a (prev <= a < next) and the
  // observations that own them
  void tooth(double a, double& prev, double& next, int& prev_t, int& next_t) const {
    prev = -1e300; next = 1e300; prev_t = -1; next_t = -1;
    for (int t = lo; t < hi; t++) {
      double N = std::exp(a + o[t]) / (1.0 - alpha);
      int K = (int)std::floor(N + 1e-9);
      if (K >= 1) { double b = std::log((double)K * (1.0 - alpha)) - o[t]; if (b > prev) { prev = b; prev_t = t; } }
      double b2 = std::log((double)(K + 1) * (1.0 - alpha)) - o[t]; if (b2 < next) { next = b2; next_t = t; }
    }
  }
  void tooth(double a, double& prev, double& next) const { int pt, nt; tooth(a, prev, next, pt, nt); }
};
}  // namespace

// window (log scale) around the peak within which breakpoints are enumerated
// (extended once on a side where the best breakpoint sits near the edge), and
// the offset below a breakpoint at which its left limit is taken
static const double BP_WINDOW = 0.15;
static const double BP_EPS    = 1e-8;
static const double GR = 0.6180339887498949;

// golden section for a maximum on [L, H]: at most `iters` iterations, stopping
// once the bracket is narrower than 1e-6 on the log scale (an intercept error
// of that size moves the objective by less than 1e-9); returns the value, a_at
// the location
static double golden_max(const Unit& U, double L, double H, int iters, double& a_at) {
  double c = H - GR * (H - L), d = L + GR * (H - L), fc = U.ll(c), fd = U.ll(d);
  for (int it = 0; it < iters && (H - L) > 1e-6; it++) {
    if (fc < fd) { L = c; c = d; fc = fd; d = L + GR * (H - L); fd = U.ll(d); }
    else         { H = d; d = c; fd = fc; c = H - GR * (H - L); fc = U.ll(c); }
  }
  a_at = 0.5 * (L + H);
  return U.ll(a_at);
}

// best left limit over the breakpoints inside (WL, WH]; best_t = its observation
static void scan_breakpoints(const Unit& U, double WL, double WH, double L0,
                             double& best_a, double& best_v, bool& best_bp, int& best_t) {
  // a breakpoint at ceiling k drops the objective by about (1-alpha)^k: only
  // ceilings whose jump exceeds 1e-7 are worth visiting
  double kcap = std::log(1e-7) / std::log1p(-U.alpha);
  int kmax_teeth = (kcap > (double)U.max_support) ? U.max_support : (int)std::floor(kcap);
  for (int t = U.lo; t < U.hi; t++) {
    double nlo = std::exp(WL + U.o[t]) / (1.0 - U.alpha), nhi = std::exp(WH + U.o[t]) / (1.0 - U.alpha);
    int klo = std::max((int)std::ceil(nlo), U.Y[t] + 1), khi = (int)std::floor(nhi);
    if (khi > kmax_teeth) khi = kmax_teeth;
    for (int k = klo; k <= khi; k++) {
      double b = std::log((double)k * (1.0 - U.alpha)) - U.o[t] - BP_EPS;
      if (b <= L0 || b <= WL) continue;
      double v = U.ll(b);
      if (v > best_v) { best_v = v; best_a = b; best_bp = true; best_t = t; }
    }
  }
}

// returns the unit's maximized log-likelihood; a_star = the intercept; edge_t =
// the row whose breakpoint the intercept sits at (-1 for an interior optimum)
static double maximize_unit(const Unit& U, double a0, double& a_star, int& edge_t, int inner_it) {
  double alo = -1e300; int ymax = 0; double sy = 0.0, se = 0.0;
  edge_t = -1;
  for (int t = U.lo; t < U.hi; t++) {
    sy += U.w[t] * (double)U.Y[t]; se += U.w[t] * std::exp(U.o[t]);
    if (U.Y[t] > 0) {
      double b = std::log(U.Y[t] * (1.0 - U.alpha)) - U.o[t];   // ceiling >= count
      if (b > alo) alo = b;
      if (U.Y[t] > ymax) ymax = U.Y[t];
    }
  }
  if (ymax == 0) { a_star = -30.0; return 0.0; }
  double apois = (sy > 0.0 && se > 0.0) ? std::log(sy / se) : alo;   // Poisson unit intercept
  double L0 = alo + 1e-9;                          // just inside the feasible set (floor tolerance 1e-9)
  bool warm = R_finite(a0), best_bp = false;
  double best_a = L0, best_v = -1e300, centre = L0; int best_t = -1;
  for (int pass = 0; pass < 2; pass++) {
    // (i) peak of the smooth part: warm start, or coarse grid + golden refinement
    if (warm) centre = std::max(a0, L0);
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
        if (centre >= H - 1e-9) { L = H; H = H + 2.0; continue; }          // peak beyond the grid: extend
        if (centre <= L + 1e-9 && L > L0 + 1e-9) { H = L; L = std::max(L0, L - 2.0); continue; }
        break;
      }
      // the peak is located to 1e-6 (golden_max stops at that bracket), finer
      // than the narrowest tooth of any ceiling the guard admits, so the tooth
      // holding the peak is identified, not guessed
      double ac; double vc = golden_max(U, std::max(L0, centre - 0.5), centre + 0.5, 80, ac);
      if (vc > bv) centre = ac;
    }
    // (ii) left limits at every breakpoint within the window. The window is the
    // base width plus the widest tooth of any observation at the centre
    // (log((K+1)/K) at ceiling K), so that units with small ceilings, whose
    // teeth are wide and tall, are searched across whole teeth.
    double wmax = 0.0;
    for (int t = U.lo; t < U.hi; t++) {
      int K = (int)std::floor(std::exp(centre + U.o[t]) / (1.0 - U.alpha) + 1e-9);
      if (K < 1) K = 1;
      double wt = std::log1p(1.0 / (double)K); if (wt > wmax) wmax = wt;
    }
    double W = std::min(1.0, BP_WINDOW + wmax);
    double WL = std::max(L0, centre - W), WH = centre + W;
    best_a = centre; best_v = U.ll(centre); best_bp = false; best_t = -1;
    scan_breakpoints(U, WL - 1e-12, WH, L0, best_a, best_v, best_bp, best_t);
    if (best_a >= WH - 0.02) { scan_breakpoints(U, WH, WH + W, L0, best_a, best_v, best_bp, best_t); WH += W; }
    else if (best_a <= WL + 0.02 && WL > L0 + 1e-9) {
      double WL2 = std::max(L0, WL - W);
      scan_breakpoints(U, WL2 - 1e-12, WL, L0, best_a, best_v, best_bp, best_t); WL = WL2;
    }
    bool at_edge = (best_a >= WH - 0.02) || (best_a <= WL + 0.02 && WL > L0 + 1e-9);
    if (warm && at_edge) { warm = false; continue; }         // peak left the warm window: cold pass
    break;
  }
  if (best_v <= -1e299) { a_star = L0; return U.ll(L0); }
  // (iii) interior critical points: golden section inside the tooth holding the
  // peak and, if the best is a breakpoint, inside the tooth ending there; never
  // below the incumbent
  double prev, next, ag, vg, prev2, next2;
  U.tooth(centre, prev, next);
  double L = std::max(L0, prev + BP_EPS), H = next - BP_EPS;
  if (H > L) { vg = golden_max(U, L, H, inner_it, ag); if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; } }
  if (best_bp) {
    U.tooth(best_a - BP_EPS, prev2, next2);
    if (std::fabs(prev2 - prev) > 1e-12) {                  // a different tooth from the one just searched
      L = std::max(L0, prev2 + BP_EPS); H = best_a;
      if (H > L) { vg = golden_max(U, L, H, inner_it, ag); if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; } }
    }
  }
  // A golden-section optimum at an end of its tooth is an edge too: at the
  // right end it is the left limit of the next breakpoint; at the left end
  // (including the feasibility floor, where an observation's ceiling equals
  // its count) it sits on the breakpoint below and tracks it as beta moves.
  // Both need the breakpoint correction in the gradient. An optimum at a
  // tooth's end may also mean the smooth peak lies in the neighbouring tooth,
  // so that tooth is searched as well, and the check repeats from there.
  for (int round = 0; round < 3 && !best_bp; round++) {
    int pt, nt; U.tooth(best_a, prev, next, pt, nt);
    bool moved = false;
    if (best_a >= next - 1e-5 && nt >= 0) {
      double aH = next - BP_EPS, vH = U.ll(aH);
      if (vH >= best_v - 1e-12) { best_v = vH; best_a = aH; }
      best_bp = true; best_t = nt;
      double p3, n3; U.tooth(next + BP_EPS, p3, n3);                 // the tooth beyond the breakpoint
      double L3 = next + BP_EPS, H3 = n3 - BP_EPS;
      if (H3 > L3) { vg = golden_max(U, L3, H3, inner_it, ag);
        if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; best_t = -1; moved = true; } }
    } else if (best_a <= prev + 1e-5 && pt >= 0) {
      double aL = std::max(L0, prev + BP_EPS), vL = U.ll(aL);
      if (vL >= best_v - 1e-12) { best_v = vL; best_a = aL; }
      best_bp = true; best_t = pt;
      double p3, n3; U.tooth(prev - BP_EPS, p3, n3);                 // the tooth ending at the breakpoint
      double L3 = std::max(L0, p3 + BP_EPS), H3 = prev - BP_EPS;
      if (H3 > L3) { vg = golden_max(U, L3, H3, inner_it, ag);
        if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; best_t = -1; moved = true; } }
    }
    if (!moved) break;
  }
  a_star = best_a; edge_t = best_bp ? best_t : -1;
  return best_v;
}

static void compute_offsets(NumericVector params, NumericMatrix X, NumericVector uoff, std::vector<double>& o) {
  int n = X.nrow(), p = X.ncol();
  for (int i = 0; i < n; i++) { double e = uoff[i]; for (int j = 0; j < p; j++) e += X(i, j) * params[j]; o[i] = e; }
}

// Concentrated negative log-likelihood at (beta, logit alpha). `awarm` carries the
// unit intercepts from the previous evaluation (length n_units; length 0 = cold
// start for every unit). Returns list(nll, a, edge): the intercepts, and for
// each unit the row whose breakpoint the intercept sits at (-1 if interior),
// which the gradient needs; on an infeasible point nll = 1e10 and `a` is
// `awarm` unchanged.
// [[Rcpp::export]]
List cpb_fe_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                    NumericVector w, IntegerVector ustart, int n_units, int max_support,
                    bool truncated, int inner_it, NumericVector awarm) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  bool warm = awarm.size() == n_units;
  IntegerVector edge(n_units, -1);
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return List::create(1e10, awarm, edge);
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector a(n_units);
  double total = 0.0, astar; int et;
  for (int u = 0; u < n_units; u++) {
    if ((u & 15) == 0) Rcpp::checkUserInterrupt();
    Unit U = { Y, o, w, ustart[u], ustart[u+1], alpha, la, l1a, max_support, truncated };
    double a0 = warm ? awarm[u] : NA_REAL;
    double ll = maximize_unit(U, a0, astar, et, inner_it);
    if (ll <= -1e299) return List::create(1e10, awarm, edge);
    a[u] = astar; edge[u] = et; total += ll;
  }
  return List::create(-total, a, edge);
}

// Gradient of the concentrated negative log-likelihood with respect to
// (beta, logit alpha) at the solved intercepts `a` (see the header comment).
// [[Rcpp::export]]
NumericVector cpb_fe_grad_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                              NumericVector w, IntegerVector ustart, int n_units, int max_support,
                              bool truncated, NumericVector a, IntegerVector edge) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || a.size() != n_units || edge.size() != n_units)
    Rcpp::stop("internal error: gradient inputs do not match the design");
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector G(p + 1);
  for (int u = 0; u < n_units; u++) {
    double Ga = 0.0, Galpha = 0.0; std::vector<double> Gb(p, 0.0);
    for (int t = ustart[u]; t < ustart[u+1]; t++) {
      double de, da;
      if (!cpb_dlogpmf(Y[t], std::exp(a[u] + o[t]), alpha, la, l1a, max_support, truncated, de, da))
        return NumericVector(p + 1, NA_REAL);
      Ga += w[t] * de; Galpha += w[t] * da;
      for (int j = 0; j < p; j++) Gb[j] += w[t] * de * X(t, j);
    }
    if (edge[u] >= 0) {                                   // the intercept tracks a breakpoint
      int ts = edge[u];
      for (int j = 0; j < p; j++) Gb[j] -= Ga * X(ts, j);
      Galpha -= Ga / (1.0 - alpha);
    }
    for (int j = 0; j < p; j++) G[j] -= Gb[j];
    G[p] -= Galpha * alpha * (1.0 - alpha);                // chain rule to logit alpha
  }
  return G;
}

// Per-observation derivatives of the CPB log pmf (columns: d/d eta, d/d alpha),
// exported for verification against numerical differences of dcpb().
// [[Rcpp::export]]
NumericMatrix cpb_dlogpmf_cpp(IntegerVector y, NumericVector lam, double alpha, int max_support, bool truncated) {
  int n = y.size();
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; i++) {
    double de, da;
    if (cpb_dlogpmf(y[i], lam[i], alpha, la, l1a, max_support, truncated, de, da)) { out(i, 0) = de; out(i, 1) = da; }
    else { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; }
  }
  return out;
}

// Unit intercepts at (beta, logit alpha) from a cold search for every unit.
// [[Rcpp::export]]
NumericVector cpb_fe_intercepts_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                                    NumericVector w, IntegerVector ustart, int n_units, int max_support,
                                    bool truncated, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector fe(n_units); double astar; int et;
  for (int u = 0; u < n_units; u++) {
    Unit U = { Y, o, w, ustart[u], ustart[u+1], alpha, la, l1a, max_support, truncated };
    maximize_unit(U, NA_REAL, astar, et, inner_it);
    fe[u] = astar;
  }
  return fe;
}
