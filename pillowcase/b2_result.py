#!/usr/bin/env python3
"""
b2_result.py -- finite q=5 candidate computation (RESEARCH_LOG sec 31).
End-to-end reproduction + assertions for the truncated polygon tables.

  B = {s_A, s_B},   two self-crossings of the blue curve
        R_t(Q_{1/3}+Q_{1/5}) at pillowcase coordinates
        s_A ~ P(0.028, 1.272)   (near the gamma=0 seam)
        s_B ~ P(3.057, 4.981)   (near the gamma=pi seam)

Within the implemented active-set search, this is the unique minimal support for
which the finite statistic n-2 rank(D_tab) is 7.  The code does not prove that B
is a bounding cochain or that D_tab is the full deformed Floer differential: it
omits branch typing, repeated insertions, arbitrary higher operations, and the
all-order Maurer--Cartan/convergence argument.

Finite-table mechanism: B activates a single immersed quadrilateral whose two
support vertices
are s_A, s_B and whose two generator-vertices are the pair {g_4,g_6} carrying one
of the two finite bigon entries; the quadrilateral cancels that entry mod 2.
The script also applies its finite monogon/distinct-input triangle obstruction
screen and, crucially, checks D_big^2=D_tab^2=0 before reporting rank statistics.
"""
import itertools
from deform import (build_geometry, bigon_matrix, triangle_contributions_P, rank_f2,
                    square_entries_f2)
from deform_full import circular_window_ok
from maurer_cartan import monogon, self_polygon
from solve_b2 import deformed, entries
from polygons import polygon_through
from earring import P_point
from bigons import _tdist


def single_quad(red, blue, gens, preims_a, preims_b, window=230):
    """The mu^3 quadrilateral matrix for one crossing pair (both preimages, both
    chain orders) -- the single entry the full build_quad computes for this pair."""
    n = len(gens)
    nb = len(blue) - 1
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            cnt = 0
            for sa in preims_a:
                for sb in preims_b:
                    if not circular_window_ok(
                            [gens[i]['kB'], gens[j]['kB'],
                             sa['kA'], sa['kB'], sb['kA'], sb['kB']], nb, window):
                        continue
                    cnt += polygon_through(red, blue, gens[i], gens[j], [sa, sb],
                                           maxspan=window)
                    cnt += polygon_through(red, blue, gens[i], gens[j], [sb, sa],
                                           maxspan=window)
            M[i][j] = cnt % 2
    return M

# target pillowcase coordinates of the two support crossings (perturbation-robust
# identification -- labels/order depend on enumeration, coordinates do not)
S_A = (0.028, 1.272)
S_B = (3.057, 4.981)


def main():
    results = []

    def check(name, ok):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    print("== P(-2,3,5): finite candidate support B={s_A,s_B} (sec 31) ==\n")
    red, blue, x = build_geometry()
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    bigons = [(i, j) for i in range(n) for j in range(n) if d[i][j]]
    raw_square = square_entries_f2(d)
    print(f"finite bigon matrix: {n} generators, entries {bigons}")
    print(f"D_big^2 audit: {'PASS' if not raw_square else 'FAIL'}; "
          f"nonzero entries {raw_square}")
    check("D_big squares to zero", not raw_square)
    print(f"finite bigon statistic h_big=n-2 rank(D_big)="
          f"{n - 2 * rank_f2(d)} (rank(D_big)={rank_f2(d)})")
    check("finite bigon gate: 9 generators, rank D_big = 2, h_big = 5",
          n == 9 and rank_f2(d) == 2 and n - 2 * rank_f2(d) == 5)

    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    crossings = [(pp, pr) for pp, pr in Pcross]
    gP = [list(P_point(g['pt'])) for g in gens]

    # locate the two support crossings by coordinate
    def find(target):
        idx = min(range(len(crossings)),
                  key=lambda k: _tdist(crossings[k][0], target))
        return idx, _tdist(crossings[idx][0], target)
    iA, dA = find(S_A)
    iB, dB = find(S_B)
    print(f"support crossings: S_A=idx{iA}@P{tuple(round(v,3) for v in crossings[iA][0])} "
          f"(d {dA:.1e}), S_B=idx{iB}@P{tuple(round(v,3) for v in crossings[iB][0])} (d {dB:.1e})")
    check("both support crossings located at the target coordinates",
          dA < 0.05 and dB < 0.05)

    Tri = {a: TriP[a] for a in range(len(TriP))}
    Qm = single_quad(red, blue, gens, crossings[iA][1], crossings[iB][1])
    Quad = {frozenset((iA, iB)): Qm}

    # the truncated table matrix for B = {S_A, S_B}
    support = (iA, iB)
    M = deformed(n, d, Tri, Quad, support)
    ent = entries(n, M)
    tab_square = square_entries_f2(M)
    print(f"D_tab^2 audit: {'PASS' if not tab_square else 'FAIL'}; "
          f"nonzero entries {tab_square}")
    check("D_tab squares to zero", not tab_square)
    print(f"truncated matrix entries: {ent}; rank(D_tab)={rank_f2(M)}; "
          f"h_tab=n-2 rank(D_tab)={n - 2 * rank_f2(M)}")
    check("D_tab has the single entry {(1,0)} (the surviving bigon entry)",
          ent == [(1, 0)])
    check("finite statistic h_tab = 7", n - 2 * rank_f2(M) == 7)

    # the support crossings carry no triangles (clean cancellation)
    check("support crossings carry no triangles",
          not any(any(r) for r in Tri[iA]) and not any(any(r) for r in Tri[iB]))

    # the canceling quadrilateral hits the bigon (4,6)-type pair
    qent = entries(n, Qm) if Qm else []
    print(f"the quadrilateral table for {{S_A,S_B}} contributes: {qent}")
    check("the quadrilateral table cancels a bigon entry",
          bool(Qm) and all(d[i][j] for (i, j) in qent))

    # Finite obstruction screen: implemented monogons and distinct-input triangles.
    nb = len(blue) - 1
    mu0 = (sum(monogon(blue, p) for p in crossings[iA][1]) % 2,
           sum(monogon(blue, p) for p in crossings[iB][1]) % 2)
    mc_viol = []
    for k in range(len(crossings)):
        if k in support:
            continue
        cnt = 0
        for pa in crossings[iA][1]:
            for pb in crossings[iB][1]:
                for pc in crossings[k][1]:
                    if not circular_window_ok(
                            [pa['kA'], pa['kB'], pb['kA'], pb['kB'],
                             pc['kA'], pc['kB']], nb, 230):
                        continue
                    for perm in itertools.permutations([pa, pb, pc]):
                        cnt += self_polygon(blue, list(perm), maxspan=230)
        if cnt % 2:
            mc_viol.append(k)
    print(f"finite obstruction screen: monogons at support = {mu0}, "
          f"distinct-input triangle violations = {mc_viol}")
    check("support passes the implemented finite obstruction screen",
          mu0 == (0, 0) and not mc_viol)

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
    if all(results):
        print("\n  finite candidate B={s_A,s_B}, with "
              "s_A~P(0.028,1.272), s_B~P(3.057,4.981)")
        print("  This output does not assert that B is a bounding cochain.")
    return all(results)


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
