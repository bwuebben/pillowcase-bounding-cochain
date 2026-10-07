#!/usr/bin/env python3
r"""Remark 5.8 and Appendix A.1(v): numerical Floer complexes for the twisted torus knots T(3, n; 2, m).

For n = 7, 8, 13 and m = 1, 2, 3, the variety of the tangle of Hedden, Herald and Kirk for T(3, n) (with the
parameters (r, s) of Section 5) is traced with the continuation code of the companion paper on T(3,n), at the two
perturbations (n eps_A, eps_B, earring) = (0.3, 0.05, 0.01) and (-0.3, -0.03, 0.02), and the image is sheared by
S_{-m}(gamma, theta) = (gamma, theta + 2m gamma) (Lemma 4.6).  The complexes are computed by the method of Remark 6.12 of
that paper (earring with delta = 0): generators, bigons and the grading formula on the arc with gr(r_+) = sigma(K),
relative gradings on each circle.  In each of the 18 cases the script checks:

  * the arc carries 1 + 2M generators and no bigon (Theorem 5.1(iii));
  * there are 2j circles, n = 6j + 1 or 6j + 2, and the homology of each circle, computed from its own complex (not
    from [HHK2, Thm. 5.4]), is (1, 1, 1, 1);
  * the total is (1 + 2j + a, 2j + b, 2j + b, 2j + a), a = floor(M/2), b = ceil(M/2) (Theorem 5.1(iii)), which equals
    (1 + floor(|sigma|/4), ceil(|sigma|/4), ceil(|sigma|/4), floor(|sigma|/4)) with sigma = sigma(T(3,n)) - 2m
    (Propositions 5.5 and 5.9).

Floating-point numerics; these computations are not used in the proofs.
"""
import sys
import time

import numpy as np

from t3n_pillowcase import Tangle, trace_variety, analyse_complex, signature, decomposition, bigon_candidates, rank_f2

CASES_N = (7, 8, 13)
CASES_M = (1, 2, 3)
PERTURBATIONS = ((0.3, 0.05, 0.01), (-0.3, -0.03, 0.02))       # (n eps_A, eps_B, earring)

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def graded_H(deg, edges):
    """F_2 homology by degree; edges (i, j) are bigons from i to j."""
    n = len(deg)
    d = np.zeros((n, n), int)
    for i, j in edges:
        d[j, i] ^= 1
    H = []
    for k in range(4):
        Ck = [i for i in range(n) if deg[i] == k]
        Ckm = [i for i in range(n) if deg[i] == (k - 1) % 4]
        Ckp = [i for i in range(n) if deg[i] == (k + 1) % 4]
        r_out = rank_f2(d[np.ix_(Ckm, Ck)]) if Ck and Ckm else 0
        r_in = rank_f2(d[np.ix_(Ck, Ckp)]) if Ck and Ckp else 0
        H.append(len(Ck) - r_out - r_in)
    return tuple(H)


def circle_complex(C, eps):
    """Relative degrees and bigons of a circle: generators in one period, targets identified modulo the period."""
    Gm, G, M = C['gens'], C['allgens'], C['M']
    per = (len(M) - 1) / 5

    def key(g):
        for j, h in enumerate(Gm):
            dpos = (h['pos'] - g['pos']) % per
            if (dpos < 2 or abs(dpos - per) < 2) and h['label'] == g['label']:
                return j
        raise KeyError(g['pos'])
    edges = [(key(b['p']), key(b['q'])) for b in bigon_candidates(M, G, eps, sources=Gm) if b['bigon']]
    return C['rel'], edges


def shear(L, m):
    L = np.array(L, float)
    L[:, 1] = L[:, 1] + 2 * m * L[:, 0]
    return L


def main():
    t0 = time.time()
    n_ok = n_cases = 0
    for n in CASES_N:
        r, s = decomposition(n)
        j = n // 6
        for (cA, eB, eps) in PERTURBATIONS:
            comps0 = trace_variety(Tangle(3, n, r, s, cA / n, eB))[0]
            for m in CASES_M:
                t1 = time.time()
                M = m if n % 3 == 1 else m + 1
                a_, b_ = M // 2, -(-M // 2)
                sig = signature(3, n) - 2 * m
                comps = [dict(C, L=shear(C['L'], m)) for C in comps0]
                R = analyse_complex(comps, eps, sig)
                A = R['arc']
                deg = [(g - sig) % 4 for g in A['gr']]
                H = list(graded_H(deg, A['bigons']))
                circ = []
                for C in R['circles']:
                    rel, edges = circle_complex(C, eps)
                    circ.append(graded_H(rel, edges))
                    H = [h + 1 for h in H]
                vec = (1 + 2 * j + a_, 2 * j + b_, 2 * j + b_, 2 * j + a_)
                ab = abs(sig)
                vec59 = (1 + ab // 4, -(-ab // 4), -(-ab // 4), ab // 4)
                ok = (len(A['gens']) == 1 + 2 * M and len(A['bigons']) == 0 and len(circ) == 2 * j
                      and all(c == (1, 1, 1, 1) for c in circ) and tuple(H) == vec == vec59)
                n_cases += 1
                n_ok += ok
                print(f'    T(3,{n};2,{m}), (n eps_A, eps_B, earring) = ({cA}, {eB}, {eps}): sigma {sig}, arc generators '
                      f'{len(A["gens"])} (1 + 2M = {1 + 2 * M}), arc bigons {len(A["bigons"])}, circles {len(circ)} '
                      f'with homology {circ}, total {tuple(H)} (Theorem 5.1(iii): {vec}) {"OK" if ok else "MISMATCH"} '
                      f'({time.time() - t1:.0f} s)', flush=True)
    check(f'all {n_cases} cases agree with Theorem 5.1(iii) (paper: 18 of 18)', n_cases == 18 and n_ok == 18)
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.0f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
