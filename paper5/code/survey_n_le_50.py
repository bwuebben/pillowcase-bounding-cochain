#!/usr/bin/env python3
r"""Remark 6.12 (last part): generator counts, bigons and homology ranks of C(L0, L1) for every
T(3, n), 4 <= n <= 50, n coprime to 3, in the decomposition of Theorem 1.1, at the four perturbations

    (n eps_A, eps_B) = (0.5, 0.03), (0.3, 0.05), (-0.5, 0.03), (0.5, -0.03),   earring eps = 0.005

(128 cases).  Checked in every case: 1 + |sigma| generators; no bigon for n = 1, 2 (mod 6) and exactly one,
from an x^- to r_+, for n = 4, 5 (mod 6); homology of rank ||Delta_K||_1; and the resolution of the double points
observed at these perturbations: split for n = 4 and joined for n = 5 (mod 6), as stated in Remark 6.12.  For
n = 5 (mod 6) all four perturbations lie in the sector containing the eps_A-axis (|lambda_c eps_A| > 2 |mu_c eps_B|); near the eps_B-axis
the resolution is split (Remark 5.8, perturbed_complexes.py).  Theorem 6.10 covers both resolutions.

Floating-point numerics; these computations are not part of the proofs.
Usage: survey_n_le_50.py [NMAX]
"""
import sys
import time

from t3n_pillowcase import (Tangle, decomposition, signature, alexander_norm, trace_variety, analyse_complex)

SETTINGS = ((0.5, 0.03), (0.3, 0.05), (-0.5, 0.03), (0.5, -0.03))
EPS = 0.005


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    t0 = time.time()
    print(f"Remark 6.12: survey of T(3,n), 4 <= n <= {nmax}, gcd(n,3) = 1, four perturbations each "
          f"(floating-point numerics)")
    rows = []
    for n in range(4, nmax + 1):
        if n % 3 == 0:
            continue
        r, s = decomposition(n)
        sig, an = signature(3, n), alexander_norm(3, n)
        for cA, eB in SETTINGS:
            t1 = time.time()
            comps, _ = trace_variety(Tangle(3, n, r, s, cA / n, eB))
            R = analyse_complex(comps, EPS, sig, gradings=False)
            A = R['arc']
            H = R['ngen'] - 2 * R['rank_d']
            cls = n % 6
            rp = next(i for i, g in enumerate(A['gens']) if g['rplus'])
            if cls in (1, 2):
                ok_b = R['nbigon'] == 0
                res = '-'
            else:
                ok_b = (R['nbigon'] == 1 and len(A['bigons']) == 1 and A['bigons'][0][1] == rp
                        and A['gens'][A['bigons'][0][0]]['label'] == 'x-')
                na = len(A['gens'])
                res = {(4, 3): 'split', (4, 7): 'joined', (5, 5): 'split', (5, 9): 'joined'}.get((cls, na), '?')
            ok_r = res == {1: '-', 2: '-', 4: 'split', 5: 'joined'}[cls]     # observed at SETTINGS (Remark 6.12)
            ok = R['ngen'] == 1 + abs(sig) and H == an and ok_b and ok_r and not R['odd']
            rows.append(ok)
            print(f"  [{'PASS' if ok else 'FAIL'}] T(3,{n:2d}) (n eA, eB)=({cA:+.1f},{eB:+.2f}): generators {R['ngen']} "
                  f"(1+|sigma| = {1 + abs(sig)}), bigons {R['nbigon']}{' x- -> r+' if R['nbigon'] == 1 and ok_b else ''}, "
                  f"rank H {H} (||Delta||_1 = {an}), resolution {res}, {len(R['circles'])} circles  [{time.time() - t1:.0f}s]",
                  flush=True)
    print(f"\nSurvey: {len(rows)} cases, {sum(rows)} agree  [{time.time() - t0:.0f}s]")
    expected = 128 if nmax == 50 else None
    if expected is not None and len(rows) != expected:
        print(f"  [FAIL] number of cases {len(rows)}, Remark 6.12: {expected}")
        return False
    if expected is not None:
        print(f"  [PASS] number of cases: computed {len(rows)}; paper 128")
    return all(rows)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
