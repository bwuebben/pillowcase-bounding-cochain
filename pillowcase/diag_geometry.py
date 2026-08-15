#!/usr/bin/env python3
"""diag_geometry.py -- dump the concrete P(-2,3,5) geometry for step (e).

Prints: the 9 generators (red x blue), the 2 bigons (which generator pairs),
blue's self-crossings (T^2 count, P-orbit count), and identifies the 4 special
connector-X self-crossings (one per resolved seam circle). This is scratch
reconnaissance for the polygon/MC computation, not a battery.
"""
import math
from tangles import TAU, PI, curve, segments, tangle_sum, fiber_circles, south_twists
from resolve import (resolve, self_crossings_T2, fold_crossings, is_closed,
                     _short_edges, _seg_int_torus)
from earring import f8, P_point, fold_pairs
from bigons import (simplify, floer_rank, intersections_detailed, arc_of,
                    is_lune, contains_vertex, CORNERS, _tdist)


def build():
    s3 = segments(curve(south_twists(3)))[0]
    s5 = segments(curve(south_twists(5)))[0]
    circ = fiber_circles(s3, s5)
    blue_polys, xinfo = resolve(tangle_sum(s3, s5), circ, eps=0.05)
    assert len(blue_polys) == 1
    blue = simplify(blue_polys[0], 3e-4)
    red = simplify(f8((2, 1), eps=0.10, phi=0.25)[0], 2e-4)
    return red, blue, blue_polys[0], xinfo, circ


if __name__ == "__main__":
    red, blue, blue_raw, xinfo, circ = build()
    print(f"blue: {len(blue)} pts (simplified); {len(blue_raw)} raw")
    print(f"red:  {len(red)} pts")

    # generators
    hits = intersections_detailed(red, blue)
    print(f"\n=== {len(hits)} generators (red x blue, C+ lift) ===")
    for i, h in enumerate(hits):
        pp = P_point(h['pt'])
        print(f"  g{i}: T2 {h['pt'][0]:.4f},{h['pt'][1]:.4f}  "
              f"P {pp[0]:.4f},{pp[1]:.4f}  redE{h['kA']} blueE{h['kB']}")

    # the 2 bigons
    print("\n=== bigons (from floer_rank verbose) ===")
    n, nl, r, hf = floer_rank(red, blue, verbose=True)
    print(f"  gens {n}, lunes {nl}, rank(d) {r}, HF {hf}")

    # blue self-crossings
    xt2 = self_crossings_T2([blue])
    orbs = fold_crossings(xt2)
    print(f"\n=== blue self-crossings: {len(xt2)} on T^2, {len(orbs)} P-orbits ===")
    # connector-X: crossings internal to each resolution's U (or D) connector pair
    print("\n=== the 4 connector-X crossings (structurally special) ===")
    for x in xinfo:
        cx = self_crossings_T2([x['U'][0], x['U'][1]])
        for p in cx:
            pp = P_point(p)
            print(f"  circle {x['circle']} seam~{x['seam']:.3f}: "
                  f"connector-X at T2 {p[0]:.4f},{p[1]:.4f}  P {pp[0]:.4f},{pp[1]:.4f}")

    print("\n=== all P-orbit representatives of blue self-crossings ===")
    for k, orb in enumerate(orbs):
        pp = P_point(orb[0])
        print(f"  s{k}: P {pp[0]:.4f},{pp[1]:.4f}  (orbit size {len(orb)})")
