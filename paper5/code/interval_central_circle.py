#!/usr/bin/env python3
r"""Lemma A.2: positivity of E_j(xi, h)/xi^3 on D_j (j = 4, 5), used in Theorem 4.10 (monotonicity of phi on
the halves of the central circle).

    D_4 = {0 <= h <= 1,   0 < xi <= pi/(2+h)},      D_5 = {0 <= h <= 1/2, 0 < xi <= pi/(2-h)},

with E_4, E_5 the functions of Lemma 4.9.  Steps, as in the proof of Lemma A.2:

 1. Exact Taylor coefficients of E_j in xi (polynomials in h with rational coefficients; exact rational
    arithmetic, see trig_series.py), compared with c3, c5, c7 as printed; the even ones vanish.
 2. The product-to-sum forms E_j = S g + q, S = sin(xi h)/h, checked exactly.
 3. The bound M9 >= |d^9 E_j / d xi^9| on [0, 1/4] x [0, h_max] from the Leibniz rule, and the remainder bound
    |R|/xi^3 <= M9 (1/4)^6 / 9!.
 4. xi <= 1/4: c3 + c5 xi^2 + c7 xi^4 - |R|/xi^3 in interval arithmetic on 64 subintervals of h.
 5. xi >= 1/4: E_j/xi^3 in interval arithmetic on 64 x 64 boxes, with sin(xi h)/h = xi sinc(xi h) and the
    monotone enclosure of sinc.
 6. Independent check stated after the lemma: splitting point xi = 0.2 and adaptive bisection of the boxes.
Interval arithmetic: mpmath.iv at 30 digits, pi enclosed.
"""
import math
import sys
import time
from fractions import Fraction as Fr

from mpmath import iv, mp, mpf

from trig_series import TrigPoly, pmul, pfrom_roots, pscale, pdiv_h, pstr

iv.dps = 30
mp.dps = 30
PI = iv.pi

PAPER = {
    4: dict(M9=9.32e7, R=0.063, small=22.06, big=0.798, small_indep=22.80),
    5: dict(M9=2.44e6, R=0.0017, small=19.18, big=0.650, small_indep=19.48),
}
HMAX = {4: Fr(1), 5: Fr(1, 2)}
results = []


def check(label, ok, computed, paper):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: computed {computed}; paper {paper}")


def lo(x):
    return mp.mpf(iv.mpf(x)._mpi_[0])


# ----------------------------------------------------------------------------- exact forms
C, S = TrigPoly.cos, TrigPoly.sin


def bracket_and_q(j):
    """the bracket multiplying sin(xi h)/h, and the remaining term, as in Lemma 4.9."""
    if j == 4:
        br = C(1, 0) * C(1, 0) * 10 + S(1, 1) * S(1, 1) * 6 - C(2, 1) * C(2, 1) * 4
        q = S(0, 1) * C(1, 0) * C(1, 0) * 4 - C(0, 1) * C(1, 0) * S(3, 2) * 2
    else:
        br = S(1, 0) * S(1, 0) * 10 + C(1, -1) * C(1, -1) * 6 - C(2, -1) * C(2, -1) * 4
        q = -(S(0, 1) * S(1, 0) * S(1, 0) * 4) - C(0, 1) * S(1, 0) * C(3, -2) * 2
    return br, q


def paper_g_q(j):
    """product-to-sum forms printed in the proof of Lemma A.2: lists of (coefficient, alpha, beta)."""
    half = Fr(1, 2)
    if j == 4:
        g = [(6, 0, 0), (5, 2, 0), (-3, 2, 2), (-2, 4, 2)]
        q = [(2, 0, 1), (1, 2, 1), (-1, 2, -1), (-half, 4, 3), (-half, 2, 1), (-half, 4, 1), (-half, 2, 3)]
    else:
        g = [(6, 0, 0), (-5, 2, 0), (3, 2, -2), (-2, 4, -2)]
        q = [(-2, 0, 1), (1, 2, 1), (-1, 2, -1), (-half, 4, -1), (half, 2, -3), (-half, 4, -3), (half, 2, -1)]
    return g, q


def as_trig(lst, kind):
    out = TrigPoly()
    for c, a, b in lst:
        out = out + (C(a, b) if kind == 'c' else S(a, b)) * Fr(c)
    return out


def paper_coeffs(j):
    if j == 4:
        base = pfrom_roots(1, [-1, -2, -3])
        c3 = pscale(base, 4)
        c5 = pscale(pmul(base, [12, 8, 5]), Fr(-2, 5))
        c7 = pscale(pmul(base, [1520, 2160, 2080, 960, 273]), Fr(1, 630))
    else:
        c3 = pscale(pmul([2, -1], [11, -2]), Fr(4, 3))
        c5 = pscale(pmul([2, -1], [-116, 152, -99, 18]), Fr(2, 15))
        c7 = pscale(pmul([-2, 1], [-4656, 11504, -14128, 9456, -3497, 502]), Fr(1, 630))
    return c3, c5, c7


def taylor_E(j, order=9):
    """Taylor coefficients of E_j in xi: (h E_j) = sin(xi h) * bracket + h * q is a trigonometric polynomial."""
    br, q = bracket_and_q(j)
    hE = S(0, 1) * br + q.scale([0, 1])
    return [pdiv_h(c) for c in hE.taylor(order)]


# ----------------------------------------------------------------------------- the bound M9
def M9_bound(j, X0):
    """M9 <= X0 |g^(9)| + sum_{k>=1} C(9,k) hmax^(k-1) |g^(9-k)| + |q^(9)|, with |f^(k)| <= sum |a_i| |w_i|^k and
    each |w_i(h)| bounded by its larger value at h = 0 and h = hmax.  Exact rational arithmetic."""
    hm = HMAX[j]
    g, q = paper_g_q(j)
    wmax = lambda a, b: max(abs(Fr(a)), abs(Fr(a) + Fr(b) * hm))
    gk = lambda k: sum(abs(Fr(c)) * wmax(a, b) ** k for c, a, b in g if (k == 0 or (a, b) != (0, 0)))
    qk = lambda k: sum(abs(Fr(c)) * wmax(a, b) ** k for c, a, b in q)
    return Fr(X0) * gk(9) + sum(math.comb(9, k) * hm ** (k - 1) * gk(9 - k) for k in range(1, 10)) + qk(9)


# ----------------------------------------------------------------------------- interval evaluation
def iv_frac(x):
    x = Fr(x)
    return iv.mpf(x.numerator) / x.denominator


def horner(p, t):
    r = iv.mpf(0)
    for c in reversed(p):
        r = r * t + iv_frac(c)
    return r


def sinc_iv(y):
    """enclosure of sin(y)/y on an interval y inside [0, pi], where sinc is decreasing."""
    def sc(t):
        t = iv.mpf(t)
        if t.a > mpf('1e-6'):
            return iv.sin(t) / t
        return 1 - t ** 2 / 6 + iv.mpf([0, 1]) * t ** 4 / 120     # 1 - y^2/6 <= sinc y <= 1 - y^2/6 + y^4/120
    return iv.mpf([sc(y.b).a, sc(y.a).b])


def E_iv(j, x, h):
    s, c = iv.sin, iv.cos
    Sx = x * sinc_iv(x * h)                                        # sin(xi h)/h
    if j == 4:
        return (Sx * (10 * c(x) ** 2 + 6 * s(x + x * h) ** 2 - 4 * c(2 * x + x * h) ** 2)
                + 4 * s(x * h) * c(x) ** 2 - 2 * c(x * h) * c(x) * s(3 * x + 2 * x * h))
    return (Sx * (10 * s(x) ** 2 + 6 * c(x - x * h) ** 2 - 4 * c(2 * x - x * h) ** 2)
            - 4 * s(x * h) * s(x) ** 2 - 2 * c(x * h) * s(x) * c(3 * x - 2 * x * h))


def xmax(j, h):
    """pi/(2 +- h), rounded up slightly so that the boxes cover D_j."""
    return (math.pi / (2 + h) if j == 4 else math.pi / (2 - h)) * (1 + 1e-12)


def small_region(j, coeffs, X0, M, nh, nx=1):
    c3, c5, c7 = coeffs
    hm = float(HMAX[j])
    rem = iv_frac(M) * iv_frac(Fr(X0)) ** 6 / math.factorial(9)
    worst = None
    for i in range(nh):
        hI = iv.mpf([hm * i / nh, hm * (i + 1) / nh])
        for k in range(nx):
            xx = iv.mpf([float(X0) * k / nx, float(X0) * (k + 1) / nx]) ** 2
            val = horner(c3, hI) + horner(c5, hI) * xx + horner(c7, hI) * xx ** 2 - rem
            worst = val.a if worst is None else min(worst, val.a)
    return worst


def big_region_grid(j, X0, nh=64, nx=64):
    hm = float(HMAX[j])
    worst, nbox, bad = None, 0, 0
    for i in range(nh):
        h1, h2 = hm * i / nh, hm * (i + 1) / nh
        xtop = max(xmax(j, h1), xmax(j, h2))
        for k in range(nx):
            x1 = X0 + (xtop - X0) * k / nx
            x2 = X0 + (xtop - X0) * (k + 1) / nx
            xI, hI = iv.mpf([x1, x2]), iv.mpf([h1, h2])
            val = E_iv(j, xI, hI) / xI ** 3
            nbox += 1
            if val.a <= 0:
                bad += 1
            worst = val.a if worst is None else min(worst, val.a)
    return worst, nbox, bad


def big_region_adaptive(j, X1, nstrip=40, maxdepth=14):
    hm = float(HMAX[j])
    stack = [(X1, max(xmax(j, hm * i / nstrip), xmax(j, hm * (i + 1) / nstrip)), hm * i / nstrip, hm * (i + 1) / nstrip, 0)
             for i in range(nstrip)]
    worst, nbox, unresolved = None, 0, 0
    while stack:
        x1, x2, h1, h2, d = stack.pop()
        if x1 > max(xmax(j, h1), xmax(j, h2)):
            continue
        xI = iv.mpf([x1, x2])
        val = E_iv(j, xI, iv.mpf([h1, h2])) / xI ** 3
        nbox += 1
        if val.a > 0:
            worst = val.a if worst is None else min(worst, val.a)
        elif d < maxdepth:
            xm, hmid = (x1 + x2) / 2, (h1 + h2) / 2
            stack += [(x1, xm, h1, hmid, d + 1), (xm, x2, h1, hmid, d + 1), (x1, xm, hmid, h2, d + 1), (xm, x2, hmid, h2, d + 1)]
        else:
            unresolved += 1
    return worst, nbox, unresolved


def main():
    t0 = time.time()
    print("Lemma A.2 (interval arithmetic, mpmath %s, %d digits, pi enclosed)" % (__import__('mpmath').__version__, iv.dps))
    for j in (4, 5):
        P = PAPER[j]
        print(f"n = {j} (mod 6):  D_{j} = {{0 <= h <= {HMAX[j]}, 0 < xi <= pi/(2{'+' if j == 4 else '-'}h)}}")
        # 1. Taylor coefficients
        T = taylor_E(j)
        c3, c5, c7 = paper_coeffs(j)
        check("even Taylor coefficients (orders 0..8) vanish identically in h",
              all(not T[k] for k in (0, 2, 4, 6, 8)) and not T[1], "yes" if all(not T[k] for k in (0, 1, 2, 4, 6, 8)) else "no",
              "E_j odd in xi, starts at xi^3")
        for name, got, want in (("c3", T[3], c3), ("c5", T[5], c5), ("c7", T[7], c7)):
            check(f"{name} equals the printed polynomial (exact)", got == want, pstr(got), "as printed")
        # 2. product-to-sum forms
        br, q = bracket_and_q(j)
        g_p, q_p = paper_g_q(j)
        check("bracket = g and remaining term = q (product-to-sum forms, exact)",
              (br - as_trig(g_p, 'c')).t == {} and (q - as_trig(q_p, 's')).t == {}, "identical", "as printed")
        # 3. M9 and the remainder
        M = M9_bound(j, Fr(1, 4))
        Rb = M * Fr(1, 4) ** 6 / math.factorial(9)
        check("M9 bound (Leibniz rule, product of maxima)", float(M) <= P['M9'], f"{float(M):.6e}", f"M9 <= {P['M9']:.3g}")
        check("|R|/xi^3 <= M9 (1/4)^6 / 9!", float(Rb) <= P['R'], f"{float(Rb):.5f}", f"<= {P['R']}")
        # 4. small xi
        low = small_region(j, (c3, c5, c7), Fr(1, 4), M, 64)
        check("xi <= 1/4, 64 subintervals of h: E_j/xi^3 >=", low >= P['small'], mp.nstr(lo(low), 6), f">= {P['small']}")
        # 5. large xi
        low, nbox, bad = big_region_grid(j, 0.25)
        check("xi >= 1/4: number of boxes", nbox == 4096, nbox, "64 x 64 = 4096")
        check("xi >= 1/4: all lower endpoints positive, minimum", bad == 0 and low >= P['big'],
              f"{bad} non-positive, min {mp.nstr(lo(low), 6)}", f"min {P['big']:.3f}")
        # 6. independent check
        M1 = M9_bound(j, Fr(1, 5))
        low_s = small_region(j, (c3, c5, c7), Fr(1, 5), M1, 80, nx=4)
        low_b, nb, unres = big_region_adaptive(j, 0.1999)   # starts below 1/5, so the two regions overlap
        check("independent check, xi <= 0.2 (80 x 4 subintervals in (h, xi)): E_j/xi^3 >=", low_s >= P['small_indep'],
              mp.nstr(lo(low_s), 6), f">= {P['small_indep']}")
        check("independent check, xi >= 0.2: adaptive boxes all positive", unres == 0 and low_b > 0,
              f"{nb} boxes, {unres} unresolved, min {mp.nstr(lo(low_b), 4)}", "positive lower bounds on all boxes")
        print("  [%.0fs]" % (time.time() - t0))
    print("\nLemma A.2: %d checks, %d passed  [%.0fs]" % (len(results), sum(results), time.time() - t0))
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
