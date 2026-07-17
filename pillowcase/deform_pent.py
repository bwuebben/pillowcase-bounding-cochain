#!/usr/bin/env python3
"""
deform_pent.py -- step (e), part 6: the mu^4 PENTAGON layer of the deformed
differential (RESEARCH_LOG sec 31).

solve_b2.py proved tri+quad insufficient: canceling the (4,6) bigon via the
{S35,S36} quad leaves (3,6),(3,4),(6,5), uncancelable below pentagon order. Here we
enumerate the 3-b-vertex polygons (immersed pentagons) g_i -> g_j through ordered
crossing TRIPLES, using the length-general polygons.polygon_through. Locality: a
pentagon is a small disk, so its 3 b-vertices and 2 generators all sit in cluster C
(gens {3,4,5,6}); triples are pruned to mutually-close near-C crossings and the two
generators to near-C. P-orbit preimages are unioned (deform convention).

Output: Pent[frozenset(a,b,c)][i][j] appended to deform_pent.json (same P-orbit
indexing as deform_full.json), consumed by solve_b2.py's cubic term.
"""
import itertools
from earring import P_point
from bigons import _tdist
from polygons import polygon_through
from deform import build_geometry, bigon_matrix, triangle_contributions_P


def build_pent(red, blue, gens, crossings, gP, cluster_gens, near_cross,
               triple_tol=0.42, gen_tol=0.80):
    """Enumerate Pent per unordered near-C crossing TRIPLE. Returns
    {frozenset(a,b,c): n x n F_2 matrix} (nonzero entries only)."""
    n = len(gens)
    Pent = {}
    triples = []
    for tri in itertools.combinations(near_cross, 3):
        pts = [crossings[t][0] for t in tri]
        if (_tdist(pts[0], pts[1]) <= triple_tol and
                _tdist(pts[1], pts[2]) <= triple_tol and
                _tdist(pts[0], pts[2]) <= triple_tol):
            triples.append(tri)
    print(f"  {len(triples)} mutually-close near-C crossing triples")
    for ti, tri in enumerate(triples):
        a, b, c = tri
        pa = [crossings[a][1], crossings[b][1], crossings[c][1]]  # preimage lists
        pp = [crossings[a][0], crossings[b][0], crossings[c][0]]
        near = [k for k in cluster_gens
                if any(_tdist(pp[m], gP[k]) < gen_tol for m in range(3))]
        M = [[0] * n for _ in range(n)]
        any_ent = False
        for i in near:
            for j in near:
                if i == j:
                    continue
                cnt = 0
                # all orderings of the 3 crossings along the blue chain, all preimages
                for perm in itertools.permutations(range(3)):
                    for pre in itertools.product(*[pa[perm[0]], pa[perm[1]], pa[perm[2]]]):
                        cnt += polygon_through(red, blue, gens[i], gens[j], list(pre))
                if cnt % 2:
                    M[i][j] = 1
                    any_ent = True
        if any_ent:
            Pent[frozenset(tri)] = M
            ent = [(i, j) for i in range(n) for j in range(n) if M[i][j]]
            print(f"    Pent{{S{a},S{b},S{c}}}: {ent}")
    return Pent


if __name__ == "__main__":
    import json
    red, blue, xinfo = build_geometry()
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    gP = [P_point(g['pt']) for g in gens]
    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    crossings = [(pp, pr) for pp, pr in Pcross]

    cluster_gens = [i for i in range(n) if abs(gP[i][1] - 4.72) < 0.3 and gP[i][0] > 3.0]
    near_cross = [a for a, (pp, pr) in enumerate(crossings)
                  if any(_tdist(pp, gP[i]) < 0.75 for i in cluster_gens)]
    print(f"cluster-C gens {cluster_gens}; near-C crossings {near_cross}")
    print("enumerating mu^4 pentagons (locality-pruned)...")
    Pent = build_pent(red, blue, gens, crossings, gP, cluster_gens, near_cross)
    print(f"pentagon triples with entries: {len(Pent)}")

    out = dict(n=n,
               Pent={",".join(map(str, sorted(k))): v for k, v in Pent.items()})
    with open("deform_pent.json", "w") as f:
        json.dump(out, f)
    print("wrote deform_pent.json")
