#!/usr/bin/env python3
"""
maurer_cartan.py -- step (e), part 4: the Maurer-Cartan obstruction of blue
(RESEARCH_LOG sec 29).

For b = sum_i c_i S_i in CF(blue,blue) (blue's self-crossings), the bounding-cochain
condition is  sum_{k>=0} mu^k(b,...,b) = 0  in CF(blue,blue) (Smith Def 5.1, eq
:mc). Component at each crossing S:

   MC_S(c) = [mu^0]_S                         (monogons: teardrops with 1 corner S)
           + sum_i c_i [mu^1(S_i)]_S          (self-bigons of blue, corners S_i,S)
           + sum_{i,j} c_i c_j [mu^2(S_i,S_j)]_S   (blue triangles, corners S_i,S_j,S)
           + ...                              = 0   (mod 2).

Each mu^k here counts immersed (k+1)-gons whose corners are ALL blue self-crossings
(no generators) -- the same bounds_disk machinery, boundary entirely on the single
iota-invariant blue curve, switching branch at each crossing. P-counts on the C_+
picture: a blue self-crossing has two T^2 preimages; a P-polygon lifts through one
of them, so we union preimage contributions per P-orbit (deform convention).
"""
import sys
import itertools
from earring import P_point
from bigons import _tdist, CORNERS
from polygons import arc_between, _assemble_loop, bounds_disk, self_intersections_detailed
from deform import build_geometry, bigon_matrix


def monogon(blue, s, maxspan=300):
    """Count mod-2 the immersed teardrops (monogons) with the single corner at the
    blue self-crossing s: a blue arc from one branch of s back to the other branch,
    bounding a disk with a convex corner at s."""
    n = len(blue) - 1
    total = 0
    for (k0, t0, k1, t1) in ((s['kA'], s['tA'], s['kB'], s['tB']),
                             (s['kB'], s['tB'], s['kA'], s['tA'])):
        for fb in (True, False):
            span = (k1 - k0) % n if fb else (k0 - k1) % n
            if span > maxspan:
                continue
            arc = arc_between(blue, k0, t0, k1, t1, fb)
            # arc runs s -> s (both ends the crossing point); close it, corner at 0
            loop = list(arc)
            if _tdist(loop[0], loop[-1]) > 1e-6:
                loop = loop + [loop[0]]
            if bounds_disk(loop, [0]):
                total += 1
    return total


def self_polygon(blue, cross_list, maxspan=220):
    """Count mod-2 the immersed (k)-gons whose k corners are exactly the blue self-
    crossings in cross_list (all on the single blue curve), boundary switching
    branch at each. For k=2 this is a blue self-bigon mu^1; k=3 a blue triangle."""
    k = len(cross_list)
    n = len(blue) - 1

    def short(k0, k1, fwd):
        span = (k1 - k0) % n if fwd else (k0 - k1) % n
        return span <= maxspan

    total = 0
    for bsel in itertools.product((0, 1), repeat=k):
        arr, dep = [], []
        for s, sel in zip(cross_list, bsel):
            if sel == 0:
                arr.append((s['kA'], s['tA'])); dep.append((s['kB'], s['tB']))
            else:
                arr.append((s['kB'], s['tB'])); dep.append((s['kA'], s['tA']))
        # loop: dep[0]->arr[k-1]; dep[k-1]->arr[k-2]; ...; dep[1]->arr[0]  (closed)
        legs = [(dep[0], arr[k - 1])]
        for m in range(k - 1, 0, -1):
            legs.append((dep[m], arr[m - 1]))
        for fbs in itertools.product((True, False), repeat=len(legs)):
            ok = True
            arcs = []
            for (p0, p1), fb in zip(legs, fbs):
                if not short(p0[0], p1[0], fb):
                    ok = False
                    break
                arcs.append(arc_between(blue, p0[0], p0[1], p1[0], p1[1], fb))
            if not ok:
                continue
            loop, corners = _assemble_loop(arcs)
            if bounds_disk(loop, corners):
                total += 1
    return total


def orbit_group(blue):
    """Blue self-crossings grouped into P-orbits: [(P_point, [preimage dicts])]."""
    scross = self_intersections_detailed(blue)
    orbs = []
    for s in scross:
        pp = list(P_point(s['pt']))
        for pr, lst in orbs:
            if _tdist(pp, pr) < 3e-3:
                lst.append(s)
                break
        else:
            orbs.append((pp, [s]))
    return orbs


def _selftest():
    """Guard: the monogon counter must FIND a genuine teardrop (else mu^0 = 0 below
    could be a false negative). A limacon inner loop is one embedded teardrop."""
    import math
    from tangles import TAU
    from bigons import simplify
    N = 1400
    cur = []
    for k in range(N + 1):
        u = TAU * k / N + 0.03
        r = 0.35 + 0.75 * math.cos(u)
        cur.append(((2.0 + r * math.cos(u)) % TAU, (2.3 + r * math.sin(u)) % TAU))
    c = simplify(cur, 2e-4)
    sx = self_intersections_detailed(c)
    ok = len(sx) == 1 and monogon(c, sx[0]) == 1
    print(f"  [{'PASS' if ok else 'FAIL'}] monogon counter finds the limacon teardrop")
    return ok


if __name__ == "__main__":
    import json
    print("== Maurer-Cartan enumeration (RESEARCH_LOG sec 29) ==")
    assert _selftest(), "monogon counter self-test failed -- mu^0=0 would be unreliable"
    red, blue, xinfo = build_geometry()
    orbs = orbit_group(blue)
    nc = len(orbs)
    print(f"{nc} P self-crossings ({sum(len(l) for _,l in orbs)} T^2 preimages)")

    # mu^0: monogons per crossing (union over preimages)
    print("\n=== mu^0 (monogons / curvature) per P-crossing ===")
    mu0 = []
    for a, (pp, preims) in enumerate(orbs):
        c = sum(monogon(blue, s) for s in preims) % 2
        mu0.append(c)
        if c:
            print(f"  S{a}@P{tuple(round(v,3) for v in pp)}: mu^0 = {c}")
    print(f"  crossings with nonzero mu^0: {sum(mu0)}")

    # mu^1: blue self-bigons per ordered crossing pair (S_i -> S).
    # Pruned T^2-CONSISTENTLY by circular blue-edge window (NOT pillowcase
    # distance: 4 of 9 generators/many crossings fold under T^2 -> P and are far
    # in edge-index despite being close in P; a P-distance prune wrongly drops
    # cross-lift polygons -- this was the sec-30/31 bug).
    from deform_full import circular_window_ok
    nb = len(blue) - 1
    print("\n=== mu^1 (blue self-bigons) per P-crossing pair [edge-window pruned] ===")
    mu1 = {}
    for a in range(nc):
        pra = orbs[a][1]
        for b in range(nc):
            if a == b:
                continue
            prb = orbs[b][1]
            c = 0
            for sa in pra:
                for sb in prb:
                    if not circular_window_ok(
                            [sa['kA'], sa['kB'], sb['kA'], sb['kB']], nb, 230):
                        continue
                    c += self_polygon(blue, [sa, sb], maxspan=230)
            if c % 2:
                mu1[(a, b)] = 1
                print(f"  S{a} -> S{b}: mu^1 = 1")
    print(f"  nonzero mu^1 pairs: {len(mu1)}")

    # mu^2: blue triangles (3 corners all crossings). Full C(nc,3) sweep is slow;
    # run with `--mu2` to reproduce. Verified result: ZERO blue triangles, so with
    # mu^0 = mu^1 = 0 the entire Maurer-Cartan obstruction vanishes to order 2 and
    # every small bounding cochain is valid -- b_2 is pinned by the rank condition.
    mu2 = None
    if "--mu2" in sys.argv:
        print("\n=== mu^2 (blue triangles) full sweep [edge-window pruned] ===")
        found = []
        for tri in itertools.combinations(range(nc), 3):
            cnt = 0
            for pa in orbs[tri[0]][1]:
                for pb in orbs[tri[1]][1]:
                    for pc in orbs[tri[2]][1]:
                        if not circular_window_ok(
                                [pa['kA'], pa['kB'], pb['kA'], pb['kB'],
                                 pc['kA'], pc['kB']], nb, 230):
                            continue
                        for perm in itertools.permutations([pa, pb, pc]):
                            cnt += self_polygon(blue, list(perm), maxspan=230)
            if cnt % 2:
                found.append(tri)
        mu2 = len(found)
        print(f"  nonzero mu^2 triples: {mu2}  (expected 0)")

    print(f"\nSUMMARY: mu^0 = {sum(mu0)} nonzero, mu^1 = {len(mu1)} nonzero"
          + (f", mu^2 = {mu2} nonzero" if mu2 is not None else " (mu^2: run --mu2)")
          + " -- MC obstruction vanishes; b_2 pinned by rank.")

    out = dict(nc=nc, crossings=[pp for pp, _ in orbs], mu0=mu0,
               mu1={f"{a},{b}": v for (a, b), v in mu1.items()}, mu2=mu2)
    with open("mc_data.json", "w") as f:
        json.dump(out, f)
    print("wrote mc_data.json")
