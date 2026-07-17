#!/usr/bin/env python3
"""pert_check.py -- perturbation-stability of the b_2 = {S13,S18} result
(RESEARCH_LOG sec 31). Runs the full step-(e) pipeline (geometry -> generators ->
bigons -> triangles -> quads -> rank-1 search -> MC check) at a chosen perturbation
and reports the minimal rank-1 bounding cochain(s). The crossing/gen LABELS are
perturbation-dependent; the INVARIANT content to confirm is: 9 gens / rank d = 2 /
HF 5, and a UNIQUE minimal-size support giving rank(partial_b) = 1 (HF 7) via a
triangle-free (4,6)-type quad, with MC = 0."""
import sys, itertools
from deform import (build_geometry_p, bigon_matrix, triangle_contributions_P,
                    rank_f2)
from deform_full import build_quad, circular_window_ok
from maurer_cartan import orbit_group, self_polygon, monogon
from solve_b2 import deformed, entries
from earring import P_point


def run(be, re, rp, window=230):
    red, blue, x = build_geometry_p(be, re, rp)
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    gP = [list(P_point(g['pt'])) for g in gens]
    bigons = [(i, j) for i in range(n) for j in range(n) if d[i][j]]
    print(f"perturbation blue={be} red_eps={re} phi={rp}: {n} gens, "
          f"rank d={rank_f2(d)} HF={n-2*rank_f2(d)}, bigons {bigons}")
    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    crossings = [(pp, pr) for pp, pr in Pcross]
    Tri = {a: TriP[a] for a in range(len(TriP))}
    Quad = {frozenset(k): v for k, v in
            build_quad(red, blue, gens, crossings, gP, window=window).items()}
    print(f"  {sum(1 for M in TriP if any(any(r) for r in M))} tri-crossings, "
          f"{len(Quad)} quad-pairs")

    act = set(a for a in Tri if any(any(r) for r in Tri[a]))
    for k in Quad:
        act |= set(k)
    active = sorted(act)
    sols = []
    for r in range(1, 4):
        for combo in itertools.combinations(active, r):
            if rank_f2(deformed(n, d, Tri, Quad, combo)) == 1:
                sols.append(combo)
        if sols:
            break
    print(f"  minimal rank-1 supports (size {r}): {len(sols)}")
    nb = len(blue) - 1
    for s in sols:
        M = deformed(n, d, Tri, Quad, s)
        # MC check: mu^0 and mu^2(b,b) over this support
        mc_ok = all(sum(monogon(blue, pre) for pre in crossings[k][1]) % 2 == 0 for k in s)
        for kk in range(len(crossings)):
            if kk in s:
                continue
            cnt = 0
            lists = [crossings[a][1] for a in s] + [crossings[kk][1]]
            if len(s) == 2:
                for pa in lists[0]:
                    for pb in lists[1]:
                        for pc in lists[2]:
                            if not circular_window_ok(
                                [pa['kA'], pa['kB'], pb['kA'], pb['kB'],
                                 pc['kA'], pc['kB']], nb, window):
                                continue
                            for perm in itertools.permutations([pa, pb, pc]):
                                cnt += self_polygon(blue, list(perm), maxspan=window)
            if cnt % 2:
                mc_ok = False
                break
        tri_free = all(not any(any(r) for r in Tri[k]) for k in s)
        Pcs = [tuple(round(v, 3) for v in crossings[k][0]) for k in s]
        print(f"    b={s} {Pcs}: partial_b {entries(n, M)}  "
              f"tri-free={tri_free}  MC={'OK' if mc_ok else 'VIOLATED'}")


if __name__ == "__main__":
    run(0.05, 0.10, 0.25)
    print()
    run(0.07, 0.16, 0.40)
