#!/usr/bin/env python3
"""
deform_full.py -- step (e), part 3: the FULL deformed differential of P(-2,3,5)
as a function of the bounding cochain b (RESEARCH_LOG sec 28/29).

partial_b[i][j] = d[i][j]
                + sum_{S in b}       Tri[S][i][j]                 (mu^2, 1 b-vertex)
                + sum_{S<S' in b}    Quad[{S,S'}][i][j]           (mu^3, 2 b-vertices)
                + ...                                              (mod 2)

Tri and Quad are enumerated over the pillowcase self-crossings (P-orbits), each the
UNION over its T^2 preimages (deform.triangle_contributions_P convention). Locality
is the tractability lever: an immersed polygon is a small disk, so its b-vertices and
its two generators all sit in one cluster -- crossing tuples and generator pairs are
pruned by mutual pillowcase distance before calling polygons.polygon_through.

Caches everything to deform_full.json so the solver (`solve_b2.py`) can evaluate
partial_b for any candidate support by table lookup.
"""
import itertools
from earring import P_point
from bigons import _tdist
from polygons import polygon_through
from deform import build_geometry, bigon_matrix, triangle_contributions_P, rank_f2


def Pdist(a, b):
    return _tdist(a, b)


def circular_window_ok(indices, n, W):
    """True if all edge-indices fit within some circular arc of length <= W on a
    cycle of n edges. This is the T^2-consistent locality test: the blue boundary
    of a small immersed disk traverses a contiguous band of blue edges, so its
    anchor edges (generators' kB and the crossings' branch edges) must cluster in
    an index window -- NOT in pillowcase (gamma,theta) distance, since 4 of the 9
    generators fold under T^2 -> P and are far in edge-index despite being close in P."""
    xs = sorted(i % n for i in indices)
    if not xs:
        return True
    gaps = [(xs[(k + 1) % len(xs)] - xs[k]) % n for k in range(len(xs))]
    return (n - max(gaps)) <= W        # complement of the largest gap = span


def build_quad(red, blue, gens, crossings, gP, window=230):
    """Enumerate Quad per unordered P-crossing pair, T^2-consistently: a quad
    g_i->g_j through crossings s,s' can only bound a small disk if the anchor blue
    edges {g_i.kB, g_j.kB, s-edge, s'-edge} fit in a circular edge window. Returns
    {frozenset(a,b): n x n F_2 matrix} (nonzero entries only)."""
    n = len(gens)
    nc = len(crossings)
    nb = len(blue) - 1
    Quad = {}
    for a in range(nc):
        pa, preims_a = crossings[a]
        for b in range(a + 1, nc):
            pb, preims_b = crossings[b]
            M = [[0] * n for _ in range(n)]
            any_ent = False
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    cnt = 0
                    for sa in preims_a:
                        for sb in preims_b:
                            anchors = [gens[i]['kB'], gens[j]['kB'],
                                       sa['kA'], sa['kB'], sb['kA'], sb['kB']]
                            if not circular_window_ok(anchors, nb, window):
                                continue
                            cnt += polygon_through(red, blue, gens[i], gens[j],
                                                   [sa, sb], maxspan=window)
                            cnt += polygon_through(red, blue, gens[i], gens[j],
                                                   [sb, sa], maxspan=window)
                    if cnt % 2:
                        M[i][j] = 1
                        any_ent = True
            if any_ent:
                Quad[frozenset((a, b))] = M
    return Quad


if __name__ == "__main__":
    import json
    red, blue, xinfo = build_geometry()
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    gP = [list(P_point(g['pt'])) for g in gens]
    print(f"{n} generators, bigons {[(i,j) for i in range(n) for j in range(n) if d[i][j]]}")

    # triangles: the validated full sweep (P-orbit indexed), NOT distance-pruned
    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    crossings = [(pp, preims) for (pp, preims) in Pcross]
    ntri = sum(1 for M in TriP if any(any(r) for r in M))
    print(f"{len(crossings)} P self-crossings; {ntri} carry triangles")

    print("enumerating quadrilaterals (locality-pruned)...")
    Quad = build_quad(red, blue, gens, crossings, gP)
    print(f"  crossing-pairs with quadrilaterals: {len(Quad)}")
    for key in sorted(Quad, key=lambda k: sorted(k)):
        M = Quad[key]
        ent = [(i, j) for i in range(n) for j in range(n) if M[i][j]]
        a, b = sorted(key)
        print(f"    Quad{{S{a}@{tuple(round(v,3) for v in crossings[a][0])}, "
              f"S{b}@{tuple(round(v,3) for v in crossings[b][0])}}}: {ent}")

    out = dict(n=n, d=d, bigons=[(i, j) for i in range(n) for j in range(n) if d[i][j]],
               gens=gP, crossings=[pp for pp, _ in crossings],
               Tri={str(a): TriP[a] for a in range(len(TriP))},
               Quad={",".join(map(str, sorted(k))): v for k, v in Quad.items()})
    with open("deform_full.json", "w") as f:
        json.dump(out, f)
    print("wrote deform_full.json")
