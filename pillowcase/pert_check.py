#!/usr/bin/env python3
"""pert_check.py -- two-perturbation check of the finite q=5 candidate tables.

Runs the implemented geometry, finite polygon tables, target-statistic search,
and finite obstruction screen. This does not compute a bounding cochain or the
full deformed Floer differential. Crossing and generator labels may change with
the perturbation; each matrix is audited for D^2=0 before a rank-derived
statistic is reported.
"""
import sys, itertools
from deform import (build_geometry_p, bigon_matrix, triangle_contributions_P,
                    rank_f2, square_entries_f2)
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
    raw_square = square_entries_f2(d)
    print(f"perturbation blue={be} red_eps={re} phi={rp}: {n} generators, "
          f"bigon entries {bigons}")
    print(f"  D_big^2 audit: {'PASS' if not raw_square else 'FAIL'}; "
          f"nonzero entries {raw_square}")
    if raw_square:
        print("  aborting before rank statistic: D_big^2 != 0")
        return
    print(f"  h_big=n-2 rank(D_big)={n-2*rank_f2(d)} "
          f"(rank(D_big)={rank_f2(d)})")
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
    audited = [(s, deformed(n, d, Tri, Quad, s)) for s in sols]
    failures = [(s, square_entries_f2(M)) for s, M in audited
                if square_entries_f2(M)]
    square_zero = [(s, M) for s, M in audited if not square_entries_f2(M)]
    print(f"  D_tab^2 audit of {len(sols)} candidate(s): "
          f"{len(square_zero)} PASS, {len(failures)} FAIL")
    for s, sq in failures:
        print(f"    FAIL support={s}: nonzero D_tab^2 entries {sq}")
    print(f"  minimal finite target supports (size {r}): {len(sols)}; "
          f"square-zero survivors: {len(square_zero)}")
    nb = len(blue) - 1
    for s, M in square_zero:
        # Finite obstruction screen: monogons and distinct-input self-triangles.
        screen_ok = all(sum(monogon(blue, pre) for pre in crossings[k][1]) % 2 == 0 for k in s)
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
                screen_ok = False
                break
        tri_free = all(not any(any(r) for r in Tri[k]) for k in s)
        Pcs = [tuple(round(v, 3) for v in crossings[k][0]) for k in s]
        print(f"    support={s} {Pcs}: D_tab^2=0; entries {entries(n, M)}; "
              f"h_tab={n-2*rank_f2(M)}; tri-free={tri_free}; "
              f"finite-screen={'PASS' if screen_ok else 'FAIL'}")


if __name__ == "__main__":
    run(0.05, 0.10, 0.25)
    print()
    run(0.07, 0.16, 0.40)
