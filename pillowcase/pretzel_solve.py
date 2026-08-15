#!/usr/bin/env python3
"""
pretzel_solve.py -- finite polygon candidate search for P(-2,3,2k+1)
(RESEARCH_LOG sec 33). For a given k this builds the piecewise-linear pillowcase
curves, a truncated matrix (bigons + triangles + distinct-support
quadrilaterals), and a finite obstruction screen (self-bigons and distinct-input
self-triangles). It searches for supports whose rank-derived statistic equals the
rigorous instanton rank.

These supports are candidates, not proved bounding cochains. In particular the
tables omit branch typing, repeated insertions, arbitrary higher operations, and
a convergence argument. A mandatory D_tab^2 check below is therefore a sanity
test, not a substitute for the missing all-order A-infinity argument.

rank I^natural(P(-2,3,q)) = q+2 rigorously for all odd q (Theorem A; see INAT
below). The finite statistic changes in different directions in the tested cases.
"""
import argparse
import itertools
from deform import (build_pretzel, bigon_matrix, triangle_contributions_P, rank_f2,
                    square_entries_f2)
from deform_full import build_quad, circular_window_ok
from maurer_cartan import orbit_group, monogon, self_polygon
from solve_b2 import deformed, entries
from earring import P_point

# rank I^natural(P(-2,3,2k+1)) = (2k+1)+2, RIGOROUS (Theorem A, RESEARCH_LOG
# sec 37 / paper2 Thm 1.1): l = sum|Delta| = q+2 (skein_alexander.py, validated)
# = u = dim Kh_r = q+2 (Manion NYJM 24 (2018) Thm 1.1) squeezes I^natural.
INAT = {k: (2 * k + 1) + 2 for k in (2, 3, 5, 6, 8, 9)}


def mu1_matrix(blue, orbs, window=230):
    """mu^1[i][j] = # blue self-bigons S_i -> S_j (mod 2), edge-window pruned."""
    nc = len(orbs); nb = len(blue) - 1
    M = [[0] * nc for _ in range(nc)]
    for i in range(nc):
        for j in range(nc):
            if i == j:
                continue
            c = 0
            for sa in orbs[i][1]:
                for sb in orbs[j][1]:
                    if not circular_window_ok([sa['kA'], sa['kB'], sb['kA'], sb['kB']], nb, window):
                        continue
                    c += self_polygon(blue, [sa, sb], maxspan=window)
            M[i][j] = c % 2
    return M


def truncated_obstruction_ok(blue, orbs, supp, mu1, window=230):
    """Apply the implementation's finite obstruction screen.

    This checks the implemented monogon/self-bigon counts and self-triangles on
    distinct untyped support points. It is not the full Maurer--Cartan equation.
    Returns (ok, detail).
    """
    nb = len(blue) - 1
    # mu^0
    for i in supp:
        if sum(monogon(blue, p) for p in orbs[i][1]) % 2:
            return False, f"monogon screen fails at S{i}"
    # mu^1(b): component at each crossing s = sum_{i in supp} mu1[i][s]
    for s in range(len(orbs)):
        v = sum(mu1[i][s] for i in supp) % 2
        # mu^2(b,b) at s: blue triangles with two corners in supp, third = s
        for (i, j) in itertools.combinations(supp, 2):
            cnt = 0
            for pa in orbs[i][1]:
                for pb in orbs[j][1]:
                    for pc in orbs[s][1]:
                        if not circular_window_ok(
                                [pa['kA'], pa['kB'], pb['kA'], pb['kB'],
                                 pc['kA'], pc['kB']], nb, window):
                            continue
                        for perm in itertools.permutations([pa, pb, pc]):
                            cnt += self_polygon(blue, list(perm), maxspan=window)
            v ^= (cnt % 2)
        if v:
            return False, f"finite obstruction screen fails at S{s}"
    return True, "ok"


def solve(k, be=0.05, re=0.16, rp=0.40, maxsupp=3, do_quad=True):
    red, blue, xi = build_pretzel(k, be, re, rp)
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    raw_square = square_entries_f2(d)
    print(f"P(-2,3,{2*k+1}): {n} generators")
    print(f"  D_big^2 audit: {'PASS' if not raw_square else 'FAIL'}; "
          f"nonzero entries {raw_square}")
    if raw_square:
        print("  aborting: no rank-derived statistic is reported for a matrix "
              "that fails D_big^2=0")
        return []
    rk0 = rank_f2(d)
    target_rank = (n - INAT[k]) // 2
    print(f"  finite bigon statistic h_big=n-2 rank(D_big)={n-2*rk0} "
          f"(rank(D_big)={rk0}); instanton benchmark={INAT[k]}")
    bigons = [(i, j) for i in range(n) for j in range(n) if d[i][j]]
    print(f"  bigon entries: {bigons}")

    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    crossings = [(pp, pr) for pp, pr in Pcross]
    gP = [list(P_point(g['pt'])) for g in gens]
    Tri = {a: TriP[a] for a in range(len(TriP))}
    Quad = {}
    if do_quad:
        print("  enumerating quadrilaterals (this is the slow step)...")
        Quad = {frozenset(kk): v for kk, v in
                build_quad(red, blue, gens, crossings, gP, window=230).items()}
    print(f"  {sum(1 for M in TriP if any(any(r) for r in M))} tri-crossings, "
          f"{len(Quad)} quad-pairs")

    orbs = orbit_group(blue)
    mu1 = mu1_matrix(blue, orbs)
    print(f"  implemented self-bigon screen nonzero pairs: "
          f"{sum(1 for i in range(len(orbs)) for j in range(len(orbs)) if mu1[i][j])}")

    active = set(a for a in Tri if any(any(r) for r in Tri[a]))
    for kk in Quad:
        active |= set(kk)
    active = sorted(active)
    print(f"  active crossings: {len(active)}")

    sols = []
    for r in range(1, maxsupp + 1):
        for combo in itertools.combinations(active, r):
            if rank_f2(deformed(n, d, Tri, Quad, combo)) == target_rank:
                ok, _ = truncated_obstruction_ok(blue, orbs, combo, mu1)
                if ok:
                    sols.append(combo)
        if sols:
            break
        print(f"  size {r}: no support passes the finite target screen")

    audited = []
    for s in sols:
        M = deformed(n, d, Tri, Quad, s)
        audited.append((s, M, square_entries_f2(M)))
    failures = [(s, sq) for s, _, sq in audited if sq]
    square_zero = [(s, M) for s, M, sq in audited if not sq]
    print(f"  D_tab^2 audit of {len(sols)} finite-screen candidate(s): "
          f"{len(square_zero)} PASS, {len(failures)} FAIL")
    for s, sq in failures:
        print(f"    FAIL support={s}: nonzero D_tab^2 entries {sq}")
    print(f"  minimal finite-screen supports with target statistic "
          f"(size {r if sols else 'n/a'}): {len(sols)}; "
          f"square-zero survivors: {len(square_zero)}")
    for s, M in square_zero[:12]:
        Pcs = [tuple(round(v, 3) for v in crossings[i][0]) for i in s]
        print(f"    support={s} {Pcs}: D_tab^2=0; entries {entries(n, M)}; "
              f"h_tab=n-2 rank(D_tab)={n-2*rank_f2(M)}")
    return sols


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the finite polygon candidate screen and mandatory D^2 audit.")
    parser.add_argument("k", nargs="?", type=int, default=3,
                        help="pretzel parameter k in P(-2,3,2k+1) (default: 3)")
    parser.add_argument("--blue-epsilon", type=float, default=0.05)
    parser.add_argument("--red-epsilon", type=float, default=0.16)
    parser.add_argument("--red-pinch", type=float, default=0.40)
    parser.add_argument("--max-support", type=int, default=3)
    parser.add_argument(
        "--triangles-only", action="store_true",
        help="skip the distinct-support quadrilateral table (appropriate for the singleton audits)")
    args = parser.parse_args()
    solve(args.k, be=args.blue_epsilon, re=args.red_epsilon,
          rp=args.red_pinch, maxsupp=args.max_support,
          do_quad=not args.triangles_only)
