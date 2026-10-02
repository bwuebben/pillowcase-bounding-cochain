#!/usr/bin/env python3
r"""Lemma A.1: the inequality used in Theorem 4.4 (monotonicity of phi on the non-central circles).

    F0 - 0.851 |F1| >= 1.648   on   Omega_+ = {0 <= delta <= pi/2,  pi/4 - delta/2 <= psi <= 3pi/4 - delta/2},

and on Omega_- = {-pi/2 <= delta <= 0, same psi-range}, where

    F0 = 10 sin^2 psi + 6 sin^2(psi + delta) - 4 sin^2(2 psi + delta),
    F1 = 4 cos(delta) sin^2 psi - 2 sin(delta) sin(psi) cos(3 psi + 2 delta).

Proof part (rigorous): Omega_+ is parametrized by (t, delta) in [0,1] x [0, pi/2], psi = pi/4 - delta/2 + pi t/2,
subdivided into 96 x 96 = 9216 boxes; on each box F0 and F1 are enclosed in interval arithmetic (mpmath.iv,
30 digits, pi enclosed in an interval) and a lower bound of F0 - 0.851 |F1| is recorded.

Independent check (rigorous point values plus Lipschitz bounds): f = F0 - 0.851 |F1| is evaluated in interval
arithmetic at the centres of a 240 x 240 grid on each of Omega_+ and Omega_-, and the Lipschitz slack
(12 pi + 6 pi * 0.851) |dt| + (8 + 10 * 0.851) |d delta| over half a cell is subtracted.

Also printed: the constant 1/(2 sin(pi/5)) bounding |lambda| in the proof of Theorem 4.4, and two
floating-point observations stated in the paper (the derivative bounds, and positivity for every c < 2).
"""
import math
import sys
import time

import numpy as np
from mpmath import iv, mp, mpf

iv.dps = 30
mp.dps = 30
PI = iv.pi
C = iv.mpf('0.851')

# values printed in the paper
PAPER = dict(boxes=9216, min_prefix='1.6482', bound=1.648, lipschitz_bound=1.637, grid=240,
             lam_bound=0.85065, deriv_bounds=(12 * math.pi, 8.0, 6 * math.pi, 10.0))

results = []


def lo(x):
    """lower endpoint of an interval, as an mpf."""
    return mp.mpf(iv.mpf(x)._mpi_[0])


def check(label, ok, computed, paper):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: computed {computed}; paper {paper}")


def F0F1(psi, d):
    """F0, F1 for interval (or float) arguments."""
    s, c = iv.sin, iv.cos
    F0 = 10 * s(psi) ** 2 + 6 * s(psi + d) ** 2 - 4 * s(2 * psi + d) ** 2
    F1 = 4 * c(d) * s(psi) ** 2 - 2 * s(d) * s(psi) * c(3 * psi + 2 * d)
    return F0, F1


def lower_f(F0, F1):
    """rigorous lower bound of F0 - 0.851 |F1| from enclosures of F0 and F1."""
    abs_hi = max(abs(F1.a), abs(F1.b))
    return (iv.mpf(F0.a) - C * iv.mpf(abs_hi)).a


def boxes(sign, n=96):
    """minimum over the n x n boxes of a lower bound of f on Omega_sign; number of boxes with bound <= 0."""
    worst, bad = None, 0
    for i in range(n):
        t = iv.mpf([mpf(i) / n, mpf(i + 1) / n])
        for j in range(n):
            d = iv.mpf([mpf(j) / n, mpf(j + 1) / n]) * (PI / 2) * sign
            psi = PI / 4 - d / 2 + t * (PI / 2)
            low = lower_f(*F0F1(psi, d))
            if low <= 0:
                bad += 1
            worst = low if worst is None else min(worst, low)
    return worst, bad


def lipschitz(sign, n=240):
    Lt = 12 * PI + 6 * PI * C
    Ld = 8 + 10 * C
    best = None
    for i in range(n):
        t = iv.mpf(2 * i + 1) / (2 * n)
        for j in range(n):
            d = sign * (iv.mpf(2 * j + 1) / (2 * n)) * (PI / 2)
            psi = PI / 4 - d / 2 + t * (PI / 2)
            low = lower_f(*F0F1(psi, d))
            best = low if best is None else min(best, low)
    slack = (Lt / (2 * n) + Ld * (PI / 2) / (2 * n)).b          # half a cell in each direction, rounded up
    return (iv.mpf(best) - iv.mpf(slack)).a, best, slack


def float_fields(t, d):
    psi = np.pi / 4 - d / 2 + np.pi * t / 2
    F0 = 10 * np.sin(psi) ** 2 + 6 * np.sin(psi + d) ** 2 - 4 * np.sin(2 * psi + d) ** 2
    F1 = 4 * np.cos(d) * np.sin(psi) ** 2 - 2 * np.sin(d) * np.sin(psi) * np.cos(3 * psi + 2 * d)
    return F0, F1


def main():
    t0 = time.time()
    print("Lemma A.1 (interval arithmetic, mpmath %s, %d digits, pi enclosed)" % (__import__('mpmath').__version__, iv.dps))

    # the constant of Theorem 4.4: max over m >= 2 of 1/(m sin(pi/(2m+1))), attained at m = 2 (n = 7)
    vals = [1 / (m * math.sin(math.pi / (2 * m + 1))) for m in range(2, 2000)]
    lam = 1 / (2 * math.sin(math.pi / 5))
    check("max_m 1/(m sin(pi/(2m+1))) is attained at m = 2 and equals 1/(2 sin(pi/5))",
          abs(max(vals) - lam) < 1e-15 and vals[0] == max(vals) and f"{lam:.5f}" == f"{PAPER['lam_bound']:.5f}",
          f"{lam:.6f}", "0.85065...")
    check("0.851 > 1/(2 sin(pi/5))", 0.851 > lam, f"{lam:.6f}", "0.851 chosen just above it")

    print("Omega_+ : %d x %d boxes" % (96, 96))
    worst, bad = boxes(+1)
    check("number of boxes", 96 * 96 == PAPER['boxes'], 96 * 96, PAPER['boxes'])
    check("no box has a non-positive lower bound", bad == 0, f"{bad} boxes", "none is negative")
    check("smallest lower bound", worst >= PAPER['bound'] and mp.nstr(lo(worst), 20).startswith(PAPER['min_prefix']),
          mp.nstr(lo(worst), 12), "1.6482...  (F0 - 0.851|F1| >= 1.648)")
    print("Omega_- : %d x %d boxes (the paper deduces this case from the symmetry (psi,delta) -> (pi-psi,-delta))" % (96, 96))
    worst_m, bad_m = boxes(-1)
    check("Omega_-: smallest lower bound >= 1.648", bad_m == 0 and worst_m >= PAPER['bound'],
          mp.nstr(lo(worst_m), 12), ">= 1.648")
    print("  [%.0fs]" % (time.time() - t0))

    print("Independent check: Lipschitz bounds and interval values at the centres of a %d x %d grid" % (PAPER['grid'], PAPER['grid']))
    # floating-point sample of the derivative bounds obtained by differentiating term by term
    rng = np.random.default_rng(0)
    t = rng.uniform(0, 1, 400000)
    d = rng.uniform(-np.pi / 2, np.pi / 2, 400000)
    h = 1e-6
    F0p, F1p = float_fields(t + h, d); F0m, F1m = float_fields(t - h, d)
    F0q, F1q = float_fields(t, d + h); F0r, F1r = float_fields(t, d - h)
    sampled = (np.max(np.abs(F0p - F0m)) / (2 * h), np.max(np.abs(F0q - F0r)) / (2 * h),
               np.max(np.abs(F1p - F1m)) / (2 * h), np.max(np.abs(F1q - F1r)) / (2 * h))
    names = ("|d_t F0| <= 12 pi", "|d_delta F0| <= 8", "|d_t F1| <= 6 pi", "|d_delta F1| <= 10")
    for name, s, b in zip(names, sampled, PAPER['deriv_bounds']):
        check(f"sampled (floating point) {name}", s <= b, f"max {s:.3f}", f"bound {b:.3f}")
    for sign, name in ((+1, "Omega_+"), (-1, "Omega_-")):
        low, best, slack = lipschitz(sign)
        check(f"{name}: f >= min over centres - slack", low >= PAPER['lipschitz_bound'],
              f"{mp.nstr(lo(best), 7)} - {mp.nstr(lo(slack), 6)} = {mp.nstr(lo(low), 6)}", f">= {PAPER['lipschitz_bound']}")
    print("  [%.0fs]" % (time.time() - t0))

    print("Floating-point observation after Lemma A.1: F0 - c|F1| > 0 on Omega_+ for every c < 2")
    tt, dd = np.meshgrid(np.linspace(0, 1, 2001), np.linspace(0, np.pi / 2, 2001))
    F0, F1 = float_fields(tt, dd)
    R = F0 / np.maximum(np.abs(F1), 1e-300)
    k = np.unravel_index(np.argmin(R), R.shape)
    # the infimum 2 is attained, e.g. at (t, delta) = (0, 0), where F0 = 4 and F1 = 2
    check("min F0/|F1| on a 2001 x 2001 grid of Omega_+ (F0 - c|F1| > 0 for all c < 2 iff this is >= 2)",
          np.min(F0) > 0 and R[k] >= 2 - 1e-12, f"min F0 = {np.min(F0):.4f}, min F0/|F1| = {R[k]:.6f} "
          f"at (t, delta) = ({tt[k]:.3f}, {dd[k]:.3f})", "positive for every c < 2")

    print("\nLemma A.1: %d checks, %d passed  [%.0fs]" % (len(results), sum(results), time.time() - t0))
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
