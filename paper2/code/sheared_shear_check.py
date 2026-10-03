#!/usr/bin/env python3
r"""Lemma 4.6(a) and Appendix A.3(iv): the regluing by t_e^h acts on the pillowcase by the shear S_h.

At points of the traced variety W_eps of the T(3,5) tangle (perturbation (5 eps_A, eps_B) = (0.3, 0.05)), the boundary
representation (a, b, c) is computed and the twisted identification (4.3) is applied:

    a -> (ba)^h a (ba)^-h,   b -> (ba)^h b (ba)^-h,   c -> c.

The pillowcase point of the result is compared, modulo the group G, with S_h(gamma, theta) = (gamma, theta - 2h gamma)
and with the shear of the opposite sign.  The paper records agreement to within 1e-12 for h = -6, ..., 2, while the
opposite shear is off by about pi.

Floating-point numerics; Lemma 4.6 is proved by hand.
"""
import sys
import time

import numpy as np

from t3n_pillowcase import Tangle, trace_variety
from sheared_complex import PI, TAU, qmul, boundary, pillow_from

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def qpow(x, k):
    r = np.tile([1.0, 0, 0, 0], (len(x), 1))
    xi = x * np.array([1, -1, -1, -1])
    for _ in range(abs(k)):
        r = qmul(r, x if k > 0 else xi)
    return r


def pdist(P1, P2):
    """Distance in the pillowcase R^2 / G between corresponding points."""
    best = None
    for sg in (1.0, -1.0):
        d = np.mod(sg * P2 - P1 + PI, TAU) - PI
        n = np.hypot(d[:, 0], d[:, 1])
        best = n if best is None else np.minimum(best, n)
    return best


def main():
    t0 = time.time()
    cA, eB = 0.3, 0.05
    comps, _ = trace_variety(Tangle(3, 5, 2, -1, cA / 5, eB), nslice=900, dimg=3e-3)
    X = np.array(comps[0]['pts'], float)
    X = X[(X[:, 0] > 0.05) & (X[:, 0] < PI - 0.05) & (X[:, 1] > 0.05) & (X[:, 1] < PI - 0.05)][::50]
    a, b, c = boundary(X, 3, 5, 2, -1, cA / 5, eB)
    P0 = pillow_from(a, b, c)
    keep = np.sin(P0[:, 0]) > 1e-3
    a, b, c, P0 = a[keep], b[keep], c[keep], P0[keep]
    print(f'{len(P0)} points of W_eps, (5 eps_A, eps_B) = ({cA}, {eB})')
    worst, opp_min = 0.0, np.inf
    for h in (-6, -5, -4, -3, -2, -1, 1, 2):
        ba = qmul(b, a)
        Wp, Wm = qpow(ba, h), qpow(ba, -h)
        P1 = pillow_from(qmul(qmul(Wp, a), Wm), qmul(qmul(Wp, b), Wm), c)
        S = np.c_[P0[:, 0], P0[:, 1] - 2 * h * P0[:, 0]]
        Sopp = np.c_[P0[:, 0], P0[:, 1] + 2 * h * P0[:, 0]]
        d, dopp = pdist(P1, S).max(), pdist(P1, Sopp).max()
        worst, opp_min = max(worst, d), min(opp_min, dopp)
        print(f'    h = {h:+d}: max distance to S_h(P) = {d:.2e}, to the opposite shear = {dopp:.2f}')
    check(f'max distance to S_h(P) over h = -6, ..., 2: {worst:.2e} (paper: within 1e-12)', worst <= 1e-12)
    check(f'the opposite shear is off by at least {opp_min:.2f} (paper: about pi)', opp_min > 2.5)
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
