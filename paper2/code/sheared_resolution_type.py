#!/usr/bin/env python3
r"""Proposition 4.10 and Remark 4.21: which resolution of the double points c_+, c_- of W occurs.

Proposition 4.10 gives d kappa_c(0) = (-tau, -1/2) at c = (pi/2, pi/2, tau), tau = +-sqrt(3)/2, so that the resolution
is joined for |eps_B| < sqrt(3) |eps_A| and split for |eps_B| > sqrt(3) |eps_A|.  The script

  (1) evaluates d kappa_c(0) = ell . d_eps Psi at c_+ and c_- with the covector ell = (-2 tau, -1) of (V5) (against the
      equations (E1), (E2); Psi = (E2, E1)), using the derivative routine of the continuation code;
  (2) traces W_eps at |eps| = 0.03 in 24 equally spaced directions, skips the four directions on the excluded lines
      (where kappa_+ or kappa_- vanishes to first order), and compares the observed type (split if and only if W_eps has
      a closed component) with the prediction.  The paper records agreement in all 20 remaining directions, the split
      ones being 75, 90, 105, 255, 270 and 285 degrees.

The shear does not change W_eps (Lemma 4.6(a)), so the answer does not depend on q.
Floating-point numerics; these computations are not used in the proofs.
"""
import math
import sys
import time

import numpy as np

from t3n_pillowcase import Tangle, trace_variety

TAU1 = math.sqrt(3) / 2
checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def main():
    t0 = time.time()
    print('(1) first-order critical values (Proposition 4.10)')
    T0 = Tangle(3, 5, 2, -1, 0.0, 0.0)
    for tau in (TAU1, -TAU1):
        D = T0.dpsi_deps(np.array([math.pi / 2, math.pi / 2, tau]))     # rows Psi_1, Psi_2; columns d/d eps_A, d/d eps_B
        ell = np.array([-1.0, -2 * tau])                                 # (-2 tau) E1 - E2, with Psi = (E2, E1)
        dk = ell @ D
        fmt = lambda v: '(' + ', '.join(f'{float(x):+.6f}' for x in v) + ')'
        check(f'tau = {tau:+.6f}: d_eps Psi_1 = {fmt(D[0])} (paper (-3 tau, 2)), d_eps Psi_2 = {fmt(D[1])} '
              f'(paper (2, -tau)), d kappa_c(0) = {fmt(dk)} (paper (-tau, -1/2))',
              np.allclose(D[0], [-3 * tau, 2], atol=1e-9) and np.allclose(D[1], [2, -tau], atol=1e-9)
              and np.allclose(dk, [-tau, -0.5], atol=1e-9))
    print('(2) 24 directions at |eps| = 0.03')
    r = 0.03
    rows, skipped = [], []
    for adeg in range(0, 360, 15):
        a = math.radians(adeg)
        eA, eB = r * math.cos(a), r * math.sin(a)
        kp, km = -TAU1 * eA - eB / 2, TAU1 * eA - eB / 2
        if min(abs(kp), abs(km)) / r < 0.08:
            skipped.append(adeg)
            print(f'    {adeg:3d} degrees: on an excluded line (kappa_+ or kappa_- vanishes to first order), skipped')
            continue
        pred = 'joined' if kp * km < 0 else 'split'
        comps, _ = trace_variety(Tangle(3, 5, 2, -1, eA, eB), nslice=900, dimg=3e-3)
        kinds = [c['kind'] for c in comps]
        obs = 'split' if any(k != 'arc' for k in kinds) else 'joined'
        rows.append((adeg, pred, obs))
        print(f'    {adeg:3d} degrees: eps = ({eA:+.4f}, {eB:+.4f}), components {kinds}: observed {obs}, '
              f'predicted {pred}', flush=True)
    agree = sum(p == o for _, p, o in rows)
    split_dirs = [a for a, _, o in rows if o == 'split']
    check(f'excluded directions {skipped} (paper: four)', skipped == [60, 120, 240, 300])
    check(f'observed type = Proposition 4.10 in {agree} of {len(rows)} directions (paper: 20 of 20)',
          len(rows) == 20 and agree == 20)
    check(f'split directions {split_dirs} (paper: 75, 90, 105, 255, 270, 285)', split_dirs == [75, 90, 105, 255, 270, 285])
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.0f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
