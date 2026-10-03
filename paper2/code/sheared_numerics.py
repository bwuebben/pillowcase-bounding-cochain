#!/usr/bin/env python3
r"""Remark 4.21, joined sectors: numerical Floer complexes of the sheared curves for P(-2,3,q).

For ten perturbations eps = (eps_A, eps_B) in the joined sectors |eps_B| < sqrt(3) |eps_A|, with |5 eps_A| between
0.002 and 0.3 and |eps_B| between 0.0003 and 0.05, the variety W_eps of the T(3,5) tangle is traced once (the shear
does not change it, Lemma 4.6(a)).  For each q in {7, 9, ..., 25, 29, 33} and the two earring parameters 0.002 and
0.0005 (240 cases), the complex of the sheared arc is computed as described in sheared_complex.py and compared with
Section 4:

  * the graded homology equals (1 + N_3, N_1, N_1, N_3) (Theorem 4.18(3)) in all 240 cases;
  * every bigon runs from a down-crossing to an up-crossing on one strand, and the number of bigons is twice the rank
    of the differential (Lemma 4.16, Proposition 4.17);
  * the degrees agree with Proposition 4.15;
  * the generator count q + 2 + 4 N_q (Proposition 4.13) is attained in 178 cases; in the other 62, which occur only at
    the larger perturbations, a line Delta_k is crossed once instead of three times (Remark 4.12).

As consistency checks the script also verifies, in every case, that the arc ends at (pi, (q+3) pi) (Proposition
4.11(1)), that its lift is embedded, that each line is crossed once (up) or three times (up, down, up), that r_+ lies
on S^0_+, and that the independently computed lift agrees with the one returned by the continuation code.

Floating-point numerics; these computations are not used in the proofs.
Usage: sheared_numerics.py [--quick]   (--quick: three perturbations, one earring, q in {7, 9, 13, 21})
"""
import math
import sys
import time

import numpy as np

from t3n_pillowcase import Tangle, trace_variety
from sheared_complex import (PI, pillow, lift, normalize_arc, shear, generators, mu_passages, sub_polyline,
                             self_crossings, bigon_test, graded_homology, predicted, Nq)

# (5 eps_A, eps_B, image step of the continuation)
PERTURBATIONS = [(+0.002, +0.0003, 2e-4), (+0.006, +0.001, 3e-4), (+0.02, +0.003, 5e-4), (+0.06, +0.01, 1e-3),
                 (+0.3, +0.05, 3e-3), (+0.3, -0.05, 3e-3), (-0.02, -0.003, 5e-4), (-0.06, -0.008, 1e-3),
                 (-0.3, -0.03, 3e-3), (-0.3, +0.03, 3e-3)]
EARRINGS = (0.002, 0.0005)
QS = (7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 29, 33)

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def analyse(Lu, q, eps):
    """The complex of the arc S_{-m}(Lu) against the earring with parameter eps."""
    L = shear(Lu, q)
    out = dict(end=tuple(int(v) for v in np.round(L[-1] / PI)))
    G = generators(L, eps)
    ir = min(range(len(G)), key=lambda j: np.hypot(*G[j]['pt']))
    G[ir]['rplus'] = True
    out['rplus_ok'] = G[ir]['sig'] == 1 and G[ir]['k'] == 0
    out['ngen'] = len(G)
    ref = math.atan2(6 - q, 1)                       # the line field (1, 6-q) = S_{-m}^{-1} l_1, Lemma 4.6(f)
    deg, deg_prop = [], []
    for j, g in enumerate(G):
        if j == ir:
            deg.append(0)
            deg_prop.append(0)
            continue
        mu = mu_passages(sub_polyline(Lu, G[ir]['pos'], g['pos']), ref)
        base = 2 * g['k'] + mu
        deg.append(base % 4 if g['sig'] > 0 else (base - 1) % 4)
        b = 2 * g['k'] if g['dir'] > 0 else 2 * g['k'] + 1            # Proposition 4.15
        deg_prop.append(b % 4 if g['sig'] > 0 else (b - 1) % 4)
    out['deg_ok'] = deg == deg_prop
    pat = {}
    for g in G:
        if g['sig'] > 0 and not g.get('rplus'):
            pat.setdefault(g['k'], []).append('u' if g['dir'] > 0 else 'd')
    out['pattern_ok'] = all(v in (['u'], ['u', 'd', 'u']) for v in pat.values())
    edges, down_up = [], True
    for ip, p in enumerate(G):
        for iq, qq in enumerate(G):
            if ip == iq or p['sig'] != qq['sig'] or p['k'] != qq['k'] or p.get('rplus') or qq.get('rplus'):
                continue
            if bigon_test(L, G, p, qq, eps)[0]:
                edges.append((ip, iq))
                down_up &= p['dir'] < 0 and qq['dir'] > 0
    H, rk = graded_homology(deg, edges, len(G))
    out.update(H=H, rank=rk, nbigon=len(edges), down_up=down_up,
               grading_drop=all((deg[i] - deg[j]) % 4 == 1 for i, j in edges))
    return out


def main():
    quick = '--quick' in sys.argv
    perts = [PERTURBATIONS[i] for i in (0, 3, 4)] if quick else PERTURBATIONS
    earrings = EARRINGS[:1] if quick else EARRINGS
    qs = (7, 9, 13, 21) if quick else QS
    t0 = time.time()
    n_cases = n_H = n_full = n_ok_struct = 0
    short_at = {}
    for (cA, eB, dimg) in perts:
        t1 = time.time()
        T = Tangle(3, 5, 2, -1, cA / 5, eB)
        comps, _ = trace_variety(T, nslice=900, dimg=dimg)
        kinds = [c['kind'] for c in comps]
        ia = kinds.index('arc')
        joined = abs(eB) < math.sqrt(3) * abs(cA / 5)
        Lmine = normalize_arc(lift(pillow(comps[ia]['pts'], 3, 5, 2, -1, cA / 5, eB)))
        Leng = normalize_arc(np.array(comps[ia]['L'], float))
        dev = float(np.abs(Lmine - Leng).max()) if Lmine.shape == Leng.shape else float('inf')
        selfx = len(self_crossings(Lmine, cell=0.02))
        print(f'(5 eps_A, eps_B) = ({cA:+.3f}, {eB:+.4f}): components {kinds}, joined sector {joined}, '
              f'lift deviation {dev:.1e}, self-crossings of the lift {selfx} ({time.time() - t1:.1f} s)', flush=True)
        check(f'  ({cA:+.3f}, {eB:+.4f}): joined sector, a single arc, embedded lift, lifts agree',
              joined and kinds == ['arc'] and selfx == 0 and dev < 1e-8)
        for eps in earrings:
            for q in qs:
                o = analyse(Lmine, q, eps)
                n_cases += 1
                full = o['ngen'] == q + 2 + 4 * Nq(q)
                n_full += full
                if not full:
                    short_at[(cA, eB)] = short_at.get((cA, eB), 0) + 1
                hok = o['H'] == predicted(q)
                n_H += hok
                struct = (o['end'] == (1, q + 3) and o['rplus_ok'] and o['pattern_ok'] and o['deg_ok'] and o['down_up']
                          and o['grading_drop'] and o['nbigon'] == 2 * o['rank'])
                n_ok_struct += struct
                print(f'    q={q:2d} eps={eps:<6g} end {o["end"]} generators {o["ngen"]:3d}/{q + 2 + 4 * Nq(q):3d} '
                      f'bigons {o["nbigon"]:2d} rank d {o["rank"]:2d} H {o["H"]} {"OK" if hok and struct else "MISMATCH"}',
                      flush=True)
    print('summary')
    check(f'graded homology = (1 + N_3, N_1, N_1, N_3) in {n_H} of {n_cases} cases'
          + ('' if quick else ' (paper: all 240)'), n_H == n_cases and (quick or n_cases == 240))
    check(f'bigons run from a down-crossing to an up-crossing on one strand, #bigons = 2 rank d, degrees as in '
          f'Proposition 4.15, end point, crossing pattern and r_+ as stated: {n_ok_struct} of {n_cases}',
          n_ok_struct == n_cases)
    if not quick:
        check(f'generator count q + 2 + 4 N_q attained in {n_full} cases, short in {n_cases - n_full} (paper: 178 and 62)',
              n_full == 178 and n_cases - n_full == 62)
        small = [(cA, eB) for (cA, eB, _) in PERTURBATIONS if abs(cA) < 0.02]
        print('    short cases by perturbation (5 eps_A, eps_B):', {k: v for k, v in sorted(short_at.items())})
        check('short counts occur only at the larger perturbations (none for |5 eps_A| < 0.02)',
              all(short_at.get(k, 0) == 0 for k in small))
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.0f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
