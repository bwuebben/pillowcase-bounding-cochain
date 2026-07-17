#!/usr/bin/env python3
"""
b2_result.py -- THE RESULT: the bounding cochain b_2 of P(-2,3,5), computed
(RESEARCH_LOG sec 31). End-to-end reproduction + assertions.

  b_2 = s_A + s_B,   two self-crossings of the blue Lagrangian
        R_t(Q_{1/3}+Q_{1/5}) at pillowcase coordinates
        s_A ~ P(0.028, 1.272)   (near the gamma=0 seam)
        s_B ~ P(3.057, 4.981)   (near the gamma=pi seam)

This is the unique MINIMAL bounding cochain for which
  rank HF((R^natural(Q_{-1/2}^hat), 0), (R_t(Q_{1/3}+Q_{1/5}), b_2)) = 7
        = rank I^natural(P(-2,3,5)),
raising the naive (b=0) rank of 5 to Kronheimer-Mrowka's 7. Smith
(arXiv:2412.06066, Thm 5.9) proved b_2 != 0 but never computed it; this is the
first explicit value of a nonzero pillowcase bounding cochain.

Mechanism: b_2 activates a single mu^3 IMMERSED QUADRILATERAL whose two b-vertices
are s_A, s_B and whose two generator-vertices are the pair {g_4,g_6} carrying one
of the two naive bigons -- the quadrilateral cancels that bigon mod 2 (rank d:
2 -> 1). The support crossings s_A, s_B carry NO other polygons (no triangles, no
monogons), so nothing else in the differential changes. Maurer-Cartan holds
automatically: blue's self-Floer operations mu^0, mu^1, mu^2 all vanish, so every
2-crossing cochain is a valid bounding cochain and the rank condition alone pins b_2.
"""
import itertools
from deform import build_geometry, bigon_matrix, triangle_contributions_P, rank_f2
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

    print("== b_2 of P(-2,3,5): the bounding cochain, computed (sec 31) ==\n")
    red, blue, x = build_geometry()
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    bigons = [(i, j) for i in range(n) for j in range(n) if d[i][j]]
    print(f"naive complex: {n} generators, bigons {bigons}, "
          f"rank d = {rank_f2(d)}, HF = {n - 2 * rank_f2(d)}")
    check("naive gate: 9 gens, rank d = 2, HF = 5 (Smith thm:main)",
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

    # the deformed differential for b_2 = {S_A, S_B}
    b = (iA, iB)
    M = deformed(n, d, Tri, Quad, b)
    ent = entries(n, M)
    print(f"deformed differential partial_b entries: {ent}, "
          f"rank = {rank_f2(M)}, HF = {n - 2 * rank_f2(M)}")
    check("partial_b has the single entry {(1,0)} (the surviving bigon)",
          ent == [(1, 0)])
    check("rank HF(b_2) = 7 = rank I^natural(P(-2,3,5))", n - 2 * rank_f2(M) == 7)

    # the support crossings carry no triangles (clean cancellation)
    check("support crossings carry no triangles",
          not any(any(r) for r in Tri[iA]) and not any(any(r) for r in Tri[iB]))

    # the canceling quadrilateral hits the bigon (4,6)-type pair
    qent = entries(n, Qm) if Qm else []
    print(f"the mu^3 quadrilateral {{S_A,S_B}} contributes: {qent}")
    check("the quadrilateral cancels a bigon (its entry is a bigon of d)",
          bool(Qm) and all(d[i][j] for (i, j) in qent))

    # Maurer-Cartan: mu^0 = 0 at the support; mu^2(b,b) = 0 -> MC = 0 exactly
    nb = len(blue) - 1
    mu0 = (sum(monogon(blue, p) for p in crossings[iA][1]) % 2,
           sum(monogon(blue, p) for p in crossings[iB][1]) % 2)
    mc_viol = []
    for k in range(len(crossings)):
        if k in b:
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
    print(f"Maurer-Cartan: mu^0 at support = {mu0}, mu^2(b,b) violations = {mc_viol}")
    check("b_2 satisfies the Maurer-Cartan equation (mu^0 = 0, mu^2(b,b) = 0)",
          mu0 == (0, 0) and not mc_viol)

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
    if all(results):
        print("\n  b_2 = s_A + s_B  with  s_A ~ P(0.028,1.272), s_B ~ P(3.057,4.981)")
        print("  -- the first explicitly computed nonzero pillowcase bounding cochain.")
    return all(results)


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
