#!/usr/bin/env python3
r"""Appendix A.1(v) (summarized in Remark 4.21), split sectors: numerical Floer complexes of the sheared arc and remnant circle for P(-2,3,q).

For five perturbations eps = (eps_A, eps_B) with |eps_B| > sqrt(3) |eps_A|, the variety W_eps of the T(3,5) tangle is an
arc and a remnant circle (Proposition 4.10).  For q = 7, 9, ..., 21 and earring parameter 0.002 (40 cases) the script
computes, as described in sheared_complex.py:

  * the sheared arc: phi-monotone, with q - 2 generators (Proposition 4.13(3)), no bigon, and homology
    (floor(K/2), ceil(K/2) - 1, ceil(K/2) - 1, floor(K/2) - 1), K = (q+1)/2 (proof of Theorem 4.18(3));
  * the remnant circle, from five periods of its sheared lift: generators in one period, relative degrees
    (Proposition 4.15), embedded bigons between generators within 1.5 periods of each other, and its homology, which
    has rank 1 in each degree (Proposition 4.17(2)), computed from its complex rather than from [HHK2, Thm. 5.4];
  * the total, which equals (1 + N_3, N_1, N_1, N_3).

The paper records: all 40 cases agree; the circle carries 4 + 4 N_q generators in 37 cases, and in the other three
(q = 17, at the three larger perturbations) a canceling pair is missing; a lift of L_0 meets the lift of the circle up to
seven times (q = 21).

Floating-point numerics; these computations are not used in the proofs.
"""
import math
import sys
import time

import numpy as np

from t3n_pillowcase import Tangle, trace_variety
from sheared_complex import (PI, normalize_arc, shear, generators, mu_passages, sub_polyline, bigon_test,
                             graded_homology, predicted, Nq)

PERTURBATIONS = [(0.0, 0.03), (0.0078, 0.029), (-0.004, -0.02), (0.0, 0.01), (0.001, 0.004)]   # (eps_A, eps_B)
QS = tuple(range(7, 22, 2))
EARRING = 0.002
EPS_BOUND = 0.031         # the bound |eps| <= 0.031 stated in Appendix A.1(v)

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def analyse_arc(L0, q, eps):
    Lu = normalize_arc(np.array(L0, float))
    L = shear(Lu, q)
    ph = L[:, 1] - L[:, 0]
    mono = bool(np.all(np.diff(ph) > -1e-12))
    G = generators(L, eps)
    ir = min(range(len(G)), key=lambda j: np.hypot(*G[j]['pt']))
    ref = math.atan2(6 - q, 1)
    deg = []
    for j, g in enumerate(G):
        if j == ir:
            deg.append(0)
            continue
        b = 2 * g['k'] + mu_passages(sub_polyline(Lu, G[ir]['pos'], g['pos']), ref)
        deg.append(b % 4 if g['sig'] > 0 else (b - 1) % 4)
    return dict(end=tuple(int(v) for v in np.round(L[-1] / PI)), mono=mono, ngen=len(G),
                H=[deg.count(i) for i in range(4)])


def analyse_circle(Lc0, q, eps, periods=5):
    Lc0 = np.array(Lc0, float)
    disp = np.round((Lc0[-1] - Lc0[0]) / PI) * PI
    closed = np.abs(Lc0[-1] - Lc0[0]).max() > 1 and np.hypot(*(Lc0[-1] - Lc0[0] - disp)) < 1e-6
    Lc = Lc0[:-1] if closed else Lc0
    n = len(Lc)
    lo = -(periods // 2)
    Mu = np.vstack([Lc + j * disp for j in range(lo, lo + periods)] + [Lc[:1] + (lo + periods) * disp])
    M = shear(Mu, q)
    G = generators(M, eps)
    base0, base1 = (-lo) * n, (-lo + 1) * n
    per_k = int(round(disp[1] / (2 * PI)))

    def cls(g):
        j = math.floor((g['pos'] - base0) / n)
        return (g['sig'], g['k'] - j * per_k, round(g['pos'] - j * n - base0, 6))
    central = [g for g in G if base0 <= g['pos'] < base1]
    index = {cls(g): i for i, g in enumerate(central)}
    ref_g = central[0]
    ref = math.atan2(6 - q, 1)
    deg = []
    for g in central:
        mu = mu_passages(sub_polyline(Mu, ref_g['pos'], g['pos']), ref) if g is not ref_g else 0
        b = 2 * (g['k'] - ref_g['k']) + mu
        deg.append(b % 4 if g['sig'] > 0 else (b - 1) % 4)
    edges = set()
    for p in central:
        for qg in G:
            if qg is p or qg['sig'] != p['sig'] or qg['k'] != p['k'] or abs(qg['pos'] - p['pos']) > 1.5 * n:
                continue
            if bigon_test(M, G, p, qg, eps)[0]:
                edges.add((index[cls(p)], index[cls(qg)]))
    H, rk = graded_homology(deg, sorted(edges), len(central))
    mult = {}
    for g in G:
        if g['sig'] > 0:
            mult[g['k']] = mult.get(g['k'], 0) + 1
    return dict(disp=tuple(int(v) for v in np.round(disp / PI)), ngen=len(central), H=list(H), rank=rk,
                nbigon=len(edges), grading_drop=all((deg[i] - deg[j]) % 4 == 1 for i, j in edges),
                maxmult=max(mult.values()) if mult else 0)


def main():
    t0 = time.time()
    n_cases = n_agree = n_full = 0
    short, maxmult = [], (0, None)
    norms = []
    for (eA, eB) in PERTURBATIONS:
        t1 = time.time()
        norms.append(math.hypot(eA, eB))
        T = Tangle(3, 5, 2, -1, eA, eB)
        comps, _ = trace_variety(T, nslice=900, dimg=2e-3)
        kinds = [c['kind'] for c in comps]
        split = abs(eB) > math.sqrt(3) * abs(eA)
        print(f'(eps_A, eps_B) = ({eA:+.4f}, {eB:+.4f}), |eps| = {math.hypot(eA, eB):.5f}: components {kinds} '
              f'({time.time() - t1:.1f} s)', flush=True)
        check(f'  split sector |eps_B| > sqrt(3)|eps_A|, and W_eps is an arc and one circle', split and sorted(kinds) == ['arc', 'circle'])
        La = comps[kinds.index('arc')]['L']
        Lc = comps[kinds.index('circle')]['L']
        for q in QS:
            K = (q + 1) // 2
            a = analyse_arc(La, q, EARRING)
            c = analyse_circle(Lc, q, EARRING)
            pred_arc = [K // 2, -(-K // 2) - 1, -(-K // 2) - 1, K // 2 - 1]
            tot = [a['H'][i] + c['H'][i] for i in range(4)]
            ok = (a['mono'] and a['ngen'] == q - 2 and a['H'] == pred_arc and a['end'] == (1, q - 1)
                  and sorted(c['H']) == [1, 1, 1, 1] and c['grading_drop'] and tot == list(predicted(q)))
            n_cases += 1
            n_agree += ok
            full = c['ngen'] == 4 + 4 * Nq(q)
            n_full += full
            if not full:
                short.append((q, (eA, eB)))
            if c['maxmult'] > maxmult[0]:
                maxmult = (c['maxmult'], q)
            print(f'    q={q:2d}: arc end {a["end"]} monotone {a["mono"]} generators {a["ngen"]}/{q - 2} H {a["H"]} | '
                  f'circle generators {c["ngen"]}/{4 + 4 * Nq(q)} bigons {c["nbigon"]} rank d {c["rank"]} H {c["H"]} '
                  f'multiplicity on a lift of L_0 {c["maxmult"]} | total {tot} {"OK" if ok else "MISMATCH"}', flush=True)
    print('summary')
    check(f'all {n_cases} cases agree with Propositions 4.13, 4.17 and Theorem 4.18 (paper: 40 of 40)',
          n_cases == 40 and n_agree == 40)
    check(f'the circle carries 4 + 4 N_q generators in {n_full} cases (paper: 37); the others are {short} '
          f'(paper: q = 17 at the three larger perturbations)',
          n_full == 37 and sorted(short) == sorted([(17, PERTURBATIONS[i]) for i in (0, 1, 2)]))
    check(f'a lift of L_0 meets the lift of the circle at most {maxmult[0]} times, attained at q = {maxmult[1]} '
          f'(paper: seven times, q = 21)', maxmult == (7, 21))
    check(f'max |eps| over the five perturbations = {max(norms):.5f} (paper: |eps| <= {EPS_BOUND})',
          max(norms) <= EPS_BOUND)
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.0f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
