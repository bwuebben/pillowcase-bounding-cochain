#!/usr/bin/env python3
"""
pretzel_solve.py -- the bounding cochain for the family P(-2,3,2k+1)
(RESEARCH_LOG sec 33). For a given k this builds the generic pillowcase complex,
the deformed differential (triangles + quadrilaterals), and the Maurer-Cartan data
(self-bigons mu^1, self-triangles mu^2), then searches for a bounding cochain b
with rank(partial_b) = (n - I^natural)/2 (i.e. deformed HF = I^natural) that also
satisfies MC. Caches per-k JSON so the search can be re-run cheaply.

I^natural(P(-2,3,2k+1)) (documented; Smith's l=u method, = sum|Alexander coeffs|):
  k=2  P(-2,3,5)=T(3,5)   : 7
  k=3  P(-2,3,7)          : 9
  k=5  P(-2,3,11)         : 13
  k=6  P(-2,3,13)         : 15
(pattern rank I^natural = 2k+3; each is naive HF + 2, a single bigon cancellation.)
"""
import sys, json, itertools
from deform import (build_pretzel, bigon_matrix, triangle_contributions_P, rank_f2)
from deform_full import build_quad, circular_window_ok
from maurer_cartan import orbit_group, monogon, self_polygon
from solve_b2 import deformed, entries
from earring import P_point

INAT = {2: 7, 3: 9, 5: 13, 6: 15}


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


def mc_ok(blue, orbs, supp, mu1, window=230):
    """Does b = sum_{i in supp} S_i satisfy Maurer-Cartan to computed order?
    mu^0 = 0 (checked), mu^1(b) then mu^2(b,b). Returns (ok, detail)."""
    nb = len(blue) - 1
    # mu^0
    for i in supp:
        if sum(monogon(blue, p) for p in orbs[i][1]) % 2:
            return False, f"mu^0 != 0 at S{i}"
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
            return False, f"MC fails at S{s}"
    return True, "ok"


def solve(k, be=0.05, re=0.16, rp=0.40, maxsupp=3, do_quad=True):
    red, blue, xi = build_pretzel(k, be, re, rp)
    gens, d = bigon_matrix(red, blue)
    n = len(gens); rk0 = rank_f2(d)
    target_rank = (n - INAT[k]) // 2
    print(f"P(-2,3,{2*k+1}): {n} gens, rank d={rk0}, naive HF={n-2*rk0}; "
          f"I^natural={INAT[k]} -> want rank(partial_b)={target_rank} (HF={n-2*target_rank})")
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
    print(f"  mu^1 (self-bigons) nonzero pairs: "
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
                ok, _ = mc_ok(blue, orbs, combo, mu1)
                if ok:
                    sols.append(combo)
        if sols:
            print(f"  minimal MC-valid rank-{target_rank} supports (size {r}): {len(sols)}")
            break
        print(f"  size {r}: no MC-valid solution (rank search only)")
    for s in sols[:12]:
        M = deformed(n, d, Tri, Quad, s)
        Pcs = [tuple(round(v, 3) for v in crossings[i][0]) for i in s]
        print(f"    b={s} {Pcs}: partial_b {entries(n, M)}  HF={n-2*rank_f2(M)}")
    return sols


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    solve(k)
