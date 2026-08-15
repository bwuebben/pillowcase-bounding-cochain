#!/usr/bin/env python3
"""diag_cancel.py -- can a polygon cancel a bigon?  Focused search for immersed
triangles (k=1) and quadrilaterals (k=2) contributing to the deformed-differential
entries of the TWO bigons of P(-2,3,5): (1,0) and (4,6). If such a polygon exists,
turning on the corresponding b-crossing(s) cancels that bigon mod 2 (rank 2 -> 1,
HF 5 -> 7). Scratch reconnaissance, not a battery."""
import itertools
from earring import P_point
from bigons import intersections_detailed, _tdist
from polygons import self_intersections_detailed, polygon_through
from deform import build_geometry, bigon_matrix

red, blue, xinfo = build_geometry()
gens, d = bigon_matrix(red, blue)
scross = self_intersections_detailed(blue)          # all 80 T^2 preimages
print(f"{len(gens)} generators, {len(scross)} T^2 self-crossings, "
      f"bigons {[(i,j) for i in range(9) for j in range(9) if d[i][j]]}")

targets = [(1, 0), (4, 6)]

print("\n=== triangles (one b-vertex) hitting a bigon entry ===")
found_tri = False
for (i, j) in targets + [(j, i) for (i, j) in targets]:
    for si, s in enumerate(scross):
        c = polygon_through(red, blue, gens[i], gens[j], [s])
        if c % 2:
            found_tri = True
            print(f"  TRIANGLE g{i}->g{j} thru S@P{tuple(round(v,4) for v in P_point(s['pt']))}")
if not found_tri:
    print("  none")

print("\n=== quadrilaterals (two b-vertices) hitting a bigon entry ===")
# restrict the crossing pair to those near the relevant generators (local disks)
def near_gens(s, idxs, tol=0.45):
    pp = P_point(s['pt'])
    return any(_tdist(pp, P_point(gens[k]['pt'])) < tol for k in idxs)

found_quad = False
for (i, j) in targets + [(j, i) for (i, j) in targets]:
    cand = [s for s in scross if near_gens(s, (i, j))]
    for sa, sb in itertools.permutations(cand, 2):
        c = polygon_through(red, blue, gens[i], gens[j], [sa, sb])
        if c % 2:
            found_quad = True
            print(f"  QUAD g{i}->g{j} thru "
                  f"S@P{tuple(round(v,4) for v in P_point(sa['pt']))} & "
                  f"S@P{tuple(round(v,4) for v in P_point(sb['pt']))}")
if not found_quad:
    print("  none among near-generator crossing pairs")
