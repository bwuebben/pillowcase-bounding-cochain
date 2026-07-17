#!/usr/bin/env python3
"""
polygons.py -- step (e), part 1: the generalized immersed-polygon counter
(RESEARCH_LOG sec 28).

Smith (arXiv:2412.06066, l.2925/2944) and CHKK: the deformed differential of
CF((L1,b1),(L2,b2)) counts immersed POLYGONS in the pillowcase whose vertices are
either generators (two of them, in CF(L1,L2)) or b-crossings (all the rest, self-
crossings of L1 or L2); the higher A_infty maps mu^n count immersed (n+1)-gons,
mu^1 = d. Everything mod 2, no signs.

This module generalizes bigons.is_lune (a 2-gon) to an arbitrary immersed k-gon:
given the closed boundary loop (a concatenation of arcs alternating between the two
Lagrangians / switching branches at self-crossings) and the list of corner indices,
`bounds_disk` decides whether the loop bounds a smooth immersed disk with convex
corners -- identical winding machinery (class 0 in H_1(T^2); w = 0 at the four
pillowcase corners; w >= 0 on every face; convex sector at each corner).

Validated against the unlink triangle: the (x,b,y) region that bigons.py test (17)
correctly REJECTS as a bigon is exactly the mu^2-triangle the bounding cochain b
activates (b = the earring pinch). This module finds precisely that one triangle.
"""
import math
from tangles import TAU, PI
from resolve import _short_edges, _seg_int_torus, is_closed
from bigons import (loop_edges, loop_class, winding_at, _corner_ok, _pt_seg_dist,
                    CORNERS, _tdist, edges_of, _seg_cross)


# ---------------------------------------------------------------------------
# self-intersections of a single curve, with edge indices (for b-crossings)
# ---------------------------------------------------------------------------
def self_intersections_detailed(poly):
    """All transverse self-intersections of a closed polyline; each returned as
    a dict with the two edge indices (kA<kB), the parameters on each, the point,
    and the two branch directions. Adjacent edges excluded."""
    E = edges_of(poly)                       # list of (k, a, b), b locally unwrapped
    n = len(poly) - 1
    closed = is_closed(poly)
    out = []
    for ii in range(len(E)):
        k1, a1, b1 = E[ii]
        for jj in range(ii + 1, len(E)):
            k2, a2, b2 = E[jj]
            dk = abs(k1 - k2)
            if dk <= 1 or (closed and dk == n - 1):
                continue
            r = _seg_cross(a1, b1, a2, b2)
            if r:
                t, u, pt = r
                p = (pt[0] % TAU, pt[1] % TAU)
                if any(_tdist(p, h['pt']) < 1e-6 for h in out):
                    continue
                out.append(dict(kA=k1, tA=t, kB=k2, tB=u, pt=p,
                                dirA=(b1[0] - a1[0], b1[1] - a1[1]),
                                dirB=(b2[0] - a2[0], b2[1] - a2[1])))
    return out


# ---------------------------------------------------------------------------
# the general immersed-disk predicate (corner-agnostic; generalizes is_lune)
# ---------------------------------------------------------------------------
def _corner_rays(loop, i):
    """Outgoing directions (away from the corner) along the two loop edges meeting
    at vertex i of a closed loop (loop[0] == loop[-1]). Handles the wrap at i==0."""
    n = len(loop) - 1                        # number of edges
    p = loop[i % n]
    nxt = loop[(i + 1) % n]
    prv = loop[(i - 1) % n]

    def d(a, b):
        return ((b[0] - a[0] + PI) % TAU - PI, (b[1] - a[1] + PI) % TAU - PI)
    return d(p, nxt), d(p, prv)              # (outgoing along next, outgoing along prev)


def bounds_disk(loop, corner_indices, verbose=False):
    """Does the closed loop (loop[0] == loop[-1]) bound a smooth immersed disk with
    convex corners exactly at `corner_indices`?  Winding machinery identical to
    bigons.is_lune, but for an arbitrary number of corners."""
    g, t, frac = loop_class(loop)
    if g != 0 or t != 0 or frac > 1e-3:
        return False
    E = loop_edges(loop)
    if len(E) < len(corner_indices):
        return False
    c0 = CORNERS[0]
    for c in CORNERS[1:]:
        if winding_at(c, E, c0) != 0:
            return False
    # w >= 0 on every face (both sides of every edge), clearance-adaptive
    nE = len(E)
    mids = [((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) for (a, b) in E]
    wmax = 0
    for i, (a, b) in enumerate(E):
        clear = 1.0
        for j, (a2, b2) in enumerate(E):
            if min((i - j) % nE, (j - i) % nE) <= 1:
                continue
            clear = min(clear, _pt_seg_dist(mids[i], a2, b2))
        L = math.hypot(b[0] - a[0], b[1] - a[1]) or 1e-12
        nx, ny = -(b[1] - a[1]) / L, (b[0] - a[0]) / L
        d = max(min(2e-3, 0.2 * L, 0.45 * clear), 1.5e-4)
        for s in (+1, -1):
            w = winding_at((mids[i][0] + s * d * nx, mids[i][1] + s * d * ny), E, c0)
            if w < 0:
                return False
            wmax = max(wmax, w)
    if wmax == 0:
        return False

    def clear_at(pp, iskip):
        c = 1.0
        for j, (a2, b2) in enumerate(E):
            if min((iskip - j) % nE, (j - iskip) % nE) <= 1:
                continue
            c = min(c, _pt_seg_dist(pp, a2, b2))
        return c

    for ci in corner_indices:
        u_dir, v_dir = _corner_rays(loop, ci)
        p = (loop[ci % (len(loop) - 1)][0] % TAU, loop[ci % (len(loop) - 1)][1] % TAU)
        if not _corner_ok(E, c0, p, u_dir, v_dir, clear_at(p, ci % nE)):
            if verbose:
                print(f"      corner {ci} rejected")
            return False
    return True


# ---------------------------------------------------------------------------
# arc extraction on a closed polyline between two (edge, param) points
# ---------------------------------------------------------------------------
def arc_between(poly, k0, t0, k1, t1, forward=True):
    """Sub-polyline of a CLOSED poly from point (k0,t0) to (k1,t1), walking in the
    direction of increasing (forward) or decreasing edge index. Points are locally
    consistent (short deltas); returns list of (gamma,theta) on T^2."""
    E = _short_edges(poly)
    n = len(poly) - 1

    def pt_on(k, tt):
        _, a, b = E[k]
        return (a[0] + tt * (b[0] - a[0]), a[1] + tt * (b[1] - a[1]))

    pts = [pt_on(k0, t0)]
    if forward:
        k = k0
        # if same edge and t1 ahead, single segment
        if k0 == k1 and t1 > t0:
            return [pt_on(k0, t0), pt_on(k1, t1)]
        while True:
            k2 = (k + 1) % n
            pts.append(poly[k2])
            if k2 == k1:
                break
            k = k2
            if k == k0:
                break
        pts.append(pt_on(k1, t1))
    else:
        k = k0
        if k0 == k1 and t1 < t0:
            return [pt_on(k0, t0), pt_on(k1, t1)]
        while True:
            pts.append(poly[k])
            k = (k - 1) % n
            if k == k1:
                break
            if k == k0:
                break
        pts.append(pt_on(k1, t1))
    return pts


def polygon_through(red, blue, x, y, cross_list, maxspan=200):
    """Count mod-2 the immersed polygons contributing to the deformed differential
    entry g_x -> g_y whose interior b-vertices are the ordered blue self-crossings
    in cross_list (each a self_intersections_detailed dict). Boundary:

        red arc(x->y) + blue arc(y -> s_k) + blue arc(s_k -> s_{k-1}) + ...
                      + blue arc(s_1 -> x),

    switching blue branch at each crossing. All branch assignments and arc
    orientations are tried; short-arc pruning keeps only local immersed disks.
    Returns the number of immersed polygons found (mod-2 caller's responsibility)."""
    import itertools
    n = len(blue) - 1

    def short(k0, k1, fwd):
        span = (k1 - k0) % n if fwd else (k0 - k1) % n
        return span <= maxspan

    k = len(cross_list)
    total = 0
    # red arc x->y (2 orientations)
    for fr in (True, False):
        red_arc = arc_between(red, x['kA'], x['tA'], y['kA'], y['tA'], fr)
        # each crossing: choose which branch is the "arrival" (0) vs "departure"
        for bsel in itertools.product((0, 1), repeat=k):
            # arrival/departure half-edges per crossing
            arr, dep = [], []
            for s, sel in zip(cross_list, bsel):
                if sel == 0:
                    arr.append((s['kA'], s['tA'])); dep.append((s['kB'], s['tB']))
                else:
                    arr.append((s['kB'], s['tB'])); dep.append((s['kA'], s['tA']))
            # blue chain: y -> arr[k-1]; dep[k-1] -> arr[k-2]; ... ; dep[0] -> x
            # waypoints in order along blue: (y) -> s_{k-1} -> ... -> s_0 -> (x)
            legs = []
            legs.append(((y['kB'], y['tB']), arr[k - 1]))
            for m in range(k - 1, 0, -1):
                legs.append((dep[m], arr[m - 1]))
            legs.append((dep[0], (x['kB'], x['tB'])))
            for fbs in itertools.product((True, False), repeat=len(legs)):
                ok = True
                blue_arcs = []
                for (p0, p1), fb in zip(legs, fbs):
                    if not short(p0[0], p1[0], fb):
                        ok = False
                        break
                    blue_arcs.append(arc_between(blue, p0[0], p0[1],
                                                 p1[0], p1[1], fb))
                if not ok:
                    continue
                loop, corners = _assemble_loop([red_arc] + blue_arcs)
                if bounds_disk(loop, corners):
                    total += 1
    return total


def _assemble_loop(arcs):
    """Concatenate arcs (each a list of pts, consecutive arcs sharing an endpoint)
    into a closed loop with loop[0] == loop[-1] and record the corner indices at
    the arc junctions."""
    loop = list(arcs[0])
    corners = [0]
    for a in arcs[1:]:
        corners.append(len(loop) - 1)
        loop += list(a)[1:]
    # loop should now return to start; make it explicitly closed
    return loop, corners


# ---------------------------------------------------------------------------
# battery (sec 28)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from tangles import curve, segments, tangle_sum, fiber_circles, south_twists
    from resolve import resolve, self_crossings_T2
    from earring import f8, arch
    from bigons import simplify, intersections_detailed

    results = []

    def check(name, ok):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    print("== the generalized immersed-polygon counter (RESEARCH_LOG sec 28) ==\n")

    # ---- (E1) sanity: bounds_disk reproduces is_lune on a synthetic 2-gon -----
    print("(E1) bounds_disk on a synthetic lens (2 circles crossing twice)")
    NA = 160
    cA = [(2.0 + 0.6 * math.cos(TAU * k / NA), 2.2 + 0.6 * math.sin(TAU * k / NA))
          for k in range(NA + 1)]
    cB = [(2.8 + 0.59 * math.cos(TAU * k / NA + 0.13),
           2.217 + 0.59 * math.sin(TAU * k / NA + 0.13)) for k in range(NA + 1)]
    hits = intersections_detailed(cA, cB)
    # build the lens loop: red arc x->y + blue arc y->x
    x, y = hits[0], hits[1]
    from polygons import arc_between as _ab
    got = 0
    for fa in (True, False):
        alpha = _ab(cA, x['kA'], x['tA'], y['kA'], y['tA'], fa)
        for fb in (True, False):
            beta = _ab(cB, y['kB'], y['tB'], x['kB'], x['tB'], fb)
            loop, corners = _assemble_loop([alpha, beta])
            if bounds_disk(loop, corners):
                got += 1
    check(f"synthetic lens: >=1 disk found by bounds_disk (got {got})", got >= 1)

    # ---- (E2) THE (x,b,y) TRIANGLE on a genuine single-curve node -------------
    # Mirrors the unlink prototype (Smith Fig UnlinkB) on exactly the code path
    # used for the real target: a curve with a GENUINE self-crossing (node b) on
    # T^2, cut twice by a transverse line -> two generators x,y bracketing the
    # node -> exactly ONE immersed (x,b,y) mu^2-triangle. (The unlink's own pinch
    # is a C_+ ^ C_- crossing between the two lifts of the earring, not a single-
    # polyline self-crossing; blue in P(-2,3,5) is a single iota-invariant curve
    # whose 40 self-crossings ARE genuine nodes, so this is the faithful test.)
    print("\n(E2) genuine figure-eight node b, cut by a transverse line:"
          " exactly ONE (x,b,y) triangle")
    Ncr = 1400
    red_fe = []
    for k in range(Ncr + 1):
        u = TAU * k / Ncr + 0.037           # jitter node off the sample grid
        red_fe.append(((1.6 + 0.9 * math.sin(u)) % TAU,
                       (2.0 + 0.6 * math.sin(2 * u)) % TAU))
    ru = simplify(red_fe, 2e-4)
    bu = simplify([((2.05) % TAU, (0.2 + 3.0 * k / 300) % TAU) for k in range(301)], 2e-4)
    gens = intersections_detailed(ru, bu)
    bcross = self_intersections_detailed(ru)          # the node (b support)
    print(f"        {len(gens)} generators, {len(bcross)} red self-crossings (nodes)")
    ntri = 0
    tri_list = []
    # triangle boundary = blue arc x->y  +  red arc y->s  +  red arc s->x
    # corners: x (gen, red^blue), y (gen, red^blue), s (red self-crossing node).
    for xi in range(len(gens)):
        for yi in range(len(gens)):
            if xi == yi:
                continue
            x, y = gens[xi], gens[yi]
            for si, s in enumerate(bcross):
                for (kIn, tIn, kOut, tOut) in ((s['kA'], s['tA'], s['kB'], s['tB']),
                                               (s['kB'], s['tB'], s['kA'], s['tA'])):
                    for fb in (True, False):
                        blue_arc = arc_between(bu, x['kB'], x['tB'], y['kB'], y['tB'], fb)
                        for f1 in (True, False):
                            red1 = arc_between(ru, y['kA'], y['tA'], kIn, tIn, f1)
                            for f2 in (True, False):
                                red2 = arc_between(ru, kOut, tOut, x['kA'], x['tA'], f2)
                                loop, corners = _assemble_loop([blue_arc, red1, red2])
                                if bounds_disk(loop, corners):
                                    ntri += 1
                                    tri_list.append((xi, yi, si))
    uniq = set((a, b) for (a, b, c) in tri_list)
    print(f"        raw triangle hits: {ntri}; distinct (x,y) gen-pairs: "
          f"{len(uniq)} {sorted(uniq)}")
    check(f"exactly one (x,b,y) triangle up to orientation (got {len(uniq)} pair)",
          len(uniq) == 1 and 1 <= ntri <= 2)

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
