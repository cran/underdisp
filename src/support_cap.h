// support_cap.h -- one rule for converting a ceiling to an int.
//
// The CPB's support runs 0..floor(n) with n = lambda / (1 - alpha), and the GEC's
// 0..floor(mu / (1 - delta)). Both diverge as the dispersion parameter approaches
// its equidispersed limit, and an optimizer reaches that limit exactly: alpha is
// carried on the logit scale, and plogis() returns exactly 1 once its argument
// passes about 40, because 1 - 4e-18 is not a representable double. The division
// then yields +Inf while the destination is an int.
//
// Converting a double outside int's range is undefined behaviour, and the
// hardware's answers disagree -- x86-64 returns INT_MIN, arm64 saturates to
// INT_MAX -- so the same fit could take different branches on different
// platforms. CRAN's M1 sanitizer build reported this on 2026-09-28, at
// cpb_pmf.h:30, cpb_fe.cpp:106 (twice) and cpb_fe.cpp:164.
//
// Every caller treats a ceiling above max_support as out of range: the pmf
// declares the point infeasible, and the fixed-effects breakpoint search never
// visits a ceiling past kmax_teeth, which is itself capped at max_support.
// Capping just past that bound before the conversion is therefore exact -- no
// case a caller can distinguish is altered -- and it replaces the hardware's
// out-of-range result with one answer on every platform.
#ifndef UNDERDISP_SUPPORT_CAP_H
#define UNDERDISP_SUPPORT_CAP_H

// x is the ceiling AFTER its floor/ceil rounding. The comparisons are written so
// that NaN takes the first branch (every comparison against NaN is false), which
// the callers reject exactly as they rejected whatever the conversion used to
// produce for it.
static inline int cap_support_int(double x, int max_support) {
  if (!(x < (double)max_support + 1.0)) return max_support + 1;  // +Inf, NaN, past the cap
  if (!(x > -1.0)) return -1;                                    // -Inf (NaN already returned)
  return (int)x;
}

#endif
