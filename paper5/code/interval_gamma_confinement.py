#!/usr/bin/env python3
r"""Lemma A.3: positivity of N_j/(h xi^{p_j}) on D_j (j = 4, 5; p_4 = 3, p_5 = 4), used in Proposition 4.11
(gamma_* < gamma < pi - gamma_* on the open halves of the central circle).

    N_4 = 4(2+h) sin(xi+xi h) cos^2 xi - (2+2h)(sin(xi+xi h) + sin(3xi+xi h)) cos^2(xi h),
    N_5 = 2(2-h) cos(xi-xi h) sin^2 xi - (cos(xi-xi h) - cos(3xi-xi h)) cos^2(xi h),

on D_4 = {0 <= h <= 1, 0 < xi <= pi/(2+h)} and D_5 = {0 <= h <= 1/2, 0 < xi <= pi/(2-h)}.  Steps, as in the proof:

 1. N_j(xi, 0) = 0 identically, and N_j = 2h cos(xi) Br_4, resp. h sin(xi) Br_5 (exact identities between
    trigonometric polynomials, with sin(alpha) = h xi sinc(alpha), alpha = xi h).
 2. Exact Taylor coefficients of N_j/h in xi: those below xi^{p_j} vanish; the leading one is as printed.
 3. xi <= 0.3: the Taylor polynomial of degrees p_j, ..., p_j + 5, plus the remainder bounded through the
    enclosure of d_xi^K d_h N_j on [0, 0.3] x [0, h], K = p_j + 6, on 64 subintervals of h.
 4. xi >= 0.3: Br_j/xi^{p_j} on boxes, starting from a 64 x 64 grid and bisecting boxes whose lower endpoint
    is not positive.
Interval arithmetic: mpmath.iv at 30 digits, pi enclosed.
"""
import math
import sys
import time
from fractions import Fraction as Fr

from mpmath import iv, mp, mpf

from trig_series import TrigPoly, pfrom_roots, pmul, pscale, pdiv_h, pstr

iv.dps = 30
mp.dps = 30

P_EXP = {4: 3, 5: 4}
HMAX = {4: 1.0, 5: 0.5}
PAPER = {4: dict(small=2.48, boxes=4100), 5: dict(small=1.28, boxes=4124)}
results = []


def check(label, ok, computed, paper):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: computed {computed}; paper {paper}")


def lo(x):
    return mp.mpf(iv.mpf(x)._mpi_[0])


C, S = TrigPoly.cos, TrigPoly.sin


def N_trig(j):
    if j == 4:
        return (S(1, 1) * C(1, 0) * C(1, 0)).scale([8, 4]) - (S(1, 1) + S(3, 1)).scale([2, 2]) * C(0, 1) * C(0, 1)
    return (C(1, -1) * S(1, 0) * S(1, 0)).scale([4, -2]) - (C(1, -1) - C(3, -1)) * C(0, 1) * C(0, 1)


def hBr_trig(j):
    """h * Br_j, written with sin(alpha) = h xi sinc(alpha), alpha = xi h."""
    sa = S(0, 1)
    if j == 4:
        s2 = S(2, 1)                                              # sin(2 xi + alpha)
        return sa.scale([2, 1]) - s2.scale([0, 1]) + (sa * sa * s2).scale([2, 2])
    s2 = S(2, -1)                                                 # sin(2 xi - alpha)
    return sa.scale([2, -1]) - s2.scale([0, 1]) + (sa * sa * s2) * 2


def leading_printed(j):
    if j == 4:
        return pscale(pmul(pfrom_roots(1, [-1, -2]), [1, 3]), Fr(4, 3))
    return pscale(pmul([2, -1], [1, 2]), Fr(2, 3))


# ----------------------------------------------------------------------------- interval evaluation
def iv_frac(x):
    x = Fr(x)
    return iv.mpf(x.numerator) / x.denominator


def horner(p, t):
    r = iv.mpf(0)
    for c in reversed(p):
        r = r * t + iv_frac(c)
    return r


def trig_iv(T, xI, hI):
    """interval enclosure of a TrigPoly at (xi, h) in the box xI x hI."""
    tot = iv.mpf(0)
    for (k, a, b), c in T.items():
        arg = (iv_frac(a) + iv_frac(b) * hI) * xI
        tot += horner(c, hI) * (iv.cos(arg) if k == 'c' else iv.sin(arg))
    return tot


def sinc_iv(y):
    """enclosure of sin(y)/y on an interval y inside [0, pi], where sinc is decreasing."""
    def sc(t):
        t = iv.mpf(t)
        if t.a > iv.mpf('1e-6').a:
            return iv.sin(t) / t
        return 1 - t ** 2 / 6 + iv.mpf([0, 1]) * t ** 4 / 120
    return iv.mpf([sc(y.b).a, sc(y.a).b])


def Br_iv(j, X, H):
    Sc = sinc_iv(X * H)
    if j == 4:
        A = 2 * X + X * H
        return (2 + H) * X * Sc - iv.sin(A) + 2 * (1 + H) * H * X ** 2 * Sc ** 2 * iv.sin(A)
    A = 2 * X - X * H
    return (2 - H) * X * Sc - iv.sin(A) + 2 * H * X ** 2 * Sc ** 2 * iv.sin(A)


def xmax(j, h):
    return math.pi / (2 + h) if j == 4 else math.pi / (2 - h)


def small_region(j, coeffs, X0=0.3, nh=64):
    p = P_EXP[j]
    K = p + 6
    poly = coeffs[p:K]                                          # degrees p, ..., p + 5 of N_j/h
    # d_xi^K (T0 + xi T1) = T0^(K) + xi T1^(K) + K T1^(K-1), where d_h N_j = T0 + xi T1
    T0, T1 = N_trig(j).dh()
    D0, D1 = T0, T1
    for _ in range(K - 1):
        D0, D1 = D0.dxi(), D1.dxi()
    D1m = D1                                                    # T1^(K-1)
    D0, D1 = D0.dxi(), D1.dxi()
    hm = HMAX[j]
    xI = iv.mpf([0, X0])
    worst = None
    for i in range(nh):
        hI = iv.mpf([hm * i / nh, hm * (i + 1) / nh])
        h0I = iv.mpf([0, hm * (i + 1) / nh])
        dK = trig_iv(D0, xI, h0I) + xI * trig_iv(D1, xI, h0I) + K * trig_iv(D1m, xI, h0I)
        MK = abs(dK).b
        rem = iv.mpf([-1, 1]) * MK * iv.mpf(X0) ** (K - p) / math.factorial(K)
        val = iv.mpf(0)
        for c in reversed(poly):
            val = val * xI + horner(c, hI)
        val = val + rem
        worst = val.a if worst is None else min(worst, val.a)
    return worst


def big_region(j, X0=0.3, nh=64, nx=64, maxdepth=12):
    p, hm = P_EXP[j], HMAX[j]
    stack = []
    for i in range(nh):
        h1, h2 = hm * i / nh, hm * (i + 1) / nh
        xtop = max(xmax(j, h1), xmax(j, h2)) * (1 + 1e-12)       # rounded up so that the boxes cover D_j
        for k in range(nx):
            stack.append((X0 + (xtop - X0) * k / nx, X0 + (xtop - X0) * (k + 1) / nx, h1, h2, 0))
    worst, nbox, bad = None, 0, 0
    while stack:
        x1, x2, h1, h2, d = stack.pop()
        if x1 > max(xmax(j, h1), xmax(j, h2)) * (1 + 1e-12):
            continue
        xI = iv.mpf([x1, x2])
        val = Br_iv(j, xI, iv.mpf([h1, h2])) / xI ** p
        nbox += 1
        if val.a > 0:
            worst = val.a if worst is None else min(worst, val.a)
        elif d < maxdepth:
            xm, hmid = (x1 + x2) / 2, (h1 + h2) / 2
            stack += [(x1, xm, h1, hmid, d + 1), (xm, x2, h1, hmid, d + 1), (x1, xm, hmid, h2, d + 1), (xm, x2, hmid, h2, d + 1)]
        else:
            bad += 1
    return worst, nbox, bad


def main():
    t0 = time.time()
    print("Lemma A.3 (interval arithmetic, mpmath %s, %d digits, pi enclosed)" % (__import__('mpmath').__version__, iv.dps))
    for j in (4, 5):
        p = P_EXP[j]
        print(f"n = {j} (mod 6), p = {p}")
        N = N_trig(j)
        check("N_j(xi, 0) = 0 identically (exact)", N.at_h0() == {}, "yes" if N.at_h0() == {} else N.at_h0(), "N_j(xi,0) = 0")
        pref = C(1, 0) * 2 if j == 4 else S(1, 0)
        check(f"N_j = {'2h cos(xi)' if j == 4 else 'h sin(xi)'} Br_j (exact)", (N - pref * hBr_trig(j)).t == {},
              "identity holds", "as printed")
        coeffs = [pdiv_h(c) for c in N.taylor(p + 6)]
        check(f"Taylor coefficients of N_j/h below xi^{p} vanish", all(not coeffs[k] for k in range(p)),
              [pstr(coeffs[k]) for k in range(p)], "vanish")
        check("leading coefficient (value at xi = 0)", coeffs[p] == leading_printed(j), pstr(coeffs[p]),
              "4/3 (h+1)(h+2)(3h+1)" if j == 4 else "2/3 (2-h)(2h+1)")
        low = small_region(j, coeffs)
        check(f"xi <= 0.3, 64 subintervals of h: N_j/(h xi^{p}) >=", low >= PAPER[j]['small'], mp.nstr(lo(low), 6),
              f">= {PAPER[j]['small']}")
        low, nbox, bad = big_region(j)
        check("xi >= 0.3: adaptive bisection from a 64 x 64 grid terminates, all lower endpoints positive",
              bad == 0 and low > 0, f"{nbox} boxes, {bad} unresolved, min Br/xi^p {mp.nstr(lo(low), 4)}",
              f"{PAPER[j]['boxes']} boxes, all positive")
        check("number of boxes", nbox == PAPER[j]['boxes'], nbox, PAPER[j]['boxes'])
        print("  [%.0fs]" % (time.time() - t0))
    print("\nLemma A.3: %d checks, %d passed  [%.0fs]" % (len(results), sum(results), time.time() - t0))
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
