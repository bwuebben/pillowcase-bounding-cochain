#!/usr/bin/env python3
"""
q7_kwz.py -- exact q=7 records, the S69 straight-through witness, and the
KWZ/type-D proof certificate for the finite PL curves (Q7_HANDOFF_2026-08-12).

The KWZ arc system must use a special puncture which no tangle arc meets.  The
q=7 blue arc meets (0,0) and (pi,0), so this module uses the CHKK top-left
choice (0,pi).  In that admissible chart the blue curve cancels to a
31-generator type-D object and its undeformed wrapped pairing with red has
dimension seven.

Smoothing the two torus lifts of each named physical node S18, S25, S69, and
S74 preserves all skeleton dots.  Every smoothing switches two same-face joins;
the four-arrow symmetric difference b is an explicit type-D Maurer--Cartan
element: D(b)=0, b^2=0, and delta+b is exactly the directly encoded smoothed
object.  Their red Hom summand dimensions are respectively 4+5, 4+5, 5+4,
and 7+2, so all four give the target wrapped dimension nine.  The S18 and S25
objects are strictly isomorphic, while their lifted component-class multiset
differs from those of S69 and S74; the four presentations therefore give
exactly three homotopy-equivalence classes.  Identifying any one of these
classes with the bounding-cochain object selected by the instanton tangle
correspondence remains separate.

Run from this directory:

    python3 q7_kwz.py            # records + witness + chart + type-D encode
    python3 q7_kwz.py --witness  # witness only
    python3 q7_kwz.py --records  # orbit records only
    python3 q7_kwz.py --encode   # type-D encode only
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction

from bigons import _tdist
from deform import build_pretzel, bigon_matrix, rank_f2
from earring import P_point
from maurer_cartan import orbit_group
from polygons import _assemble_loop, arc_between, bounds_disk
from tangles import PI, TAU

# Published default perturbation (handoff sec. 2).
DEFAULT = dict(blue_eps=0.05, red_eps=0.16, red_phi=0.40)

# CHKK Section 11.5 singles out the pillowcase corner missed by tangle
# character varieties.  Translation by this two-torsion point commutes with
# the pillowcase involution and lets the existing (0,0)-centred chart encode
# the correct special puncture.
KWZ_SPECIAL = (0.0, PI)

# Physical orbits that survive the degree-one screen, identified by pillowcase
# coordinates at DEFAULT (handoff secs. 2 and 5).
TARGETS = {
    "S69": (0.041596649200926095, 5.405395334342216),
    "S18": (0.0304313, 4.5024806),
    "S25": (3.0828526, 1.7672244),
    "S74": (3.079988, 3.561056),
}

# Session record of the straight-through self-bigon (handoff sec. 4).
WITNESS_EDGES = ((6, 905), (43, 945))
WITNESS_THROUGH_PARAM = 941.640214

FAILS = []


def check(cond, msg):
    print(("  [PASS] " if cond else "  [FAIL] ") + msg)
    if not cond:
        FAILS.append(msg)
    return cond


def loc_param(k, t):
    """Continuous edge+parameter coordinate along a closed polyline."""
    return k + t


def unwrap_path(pts):
    """Lift a T^2 polyline to R^2 by short increments."""
    if not pts:
        return []
    out = [list(pts[0])]
    for p in pts[1:]:
        dg = (p[0] - out[-1][0] + PI) % TAU - PI
        dt = (p[1] - out[-1][1] + PI) % TAU - PI
        out.append([out[-1][0] + dg, out[-1][1] + dt])
    return out


def deck_of_path(pts):
    """Integer homology class of a nearly-closed T^2 path, as (m, n) in Z^2."""
    u = unwrap_path(pts)
    if len(u) < 2:
        return (0, 0)
    dg = u[-1][0] - u[0][0]
    dt = u[-1][1] - u[0][1]
    return (int(round(dg / TAU)), int(round(dt / TAU)))


def tangent_phase(dxy):
    """Argument of a short T^2 direction, in (-pi, pi]."""
    return math.atan2(dxy[1], dxy[0])


def oriented_angle(u, v):
    """Oriented angle from direction u to direction v, in (-pi, pi]."""
    return (tangent_phase(v) - tangent_phase(u) + PI) % TAU - PI


def akaho_degrees(s):
    """Complementary degrees of the two ordered branch jumps at a crossing.

    Convention checked against the published Akaho--Joyce turning formula for
    an oriented curve in an oriented surface: the jump A -> B has degree 1 iff
    the oriented angle from dirA to dirB is positive.  The two orderings are
    complementary.  This is the convention the handoff used for S69.
    """
    ang = oriented_angle(s["dirA"], s["dirB"])
    deg_AB = 1 if ang > 0 else 0
    return {"A->B": deg_AB, "B->A": 1 - deg_AB, "angle": ang}


def branch_loop(blue, s, start_branch):
    """Closed walk that leaves s on `start_branch` and returns on the other.

    start_branch is 'A' or 'B'.  The walk follows increasing edge index from
    the departure half-edge to the arrival half-edge.
    """
    if start_branch == "A":
        k0, t0, k1, t1 = s["kA"], s["tA"], s["kB"], s["tB"]
    else:
        k0, t0, k1, t1 = s["kB"], s["tB"], s["kA"], s["tA"]
    return arc_between(blue, k0, t0, k1, t1, forward=True)


def short_deck(loop_fwd, loop_rev):
    """The shorter of the two branch-to-branch classes, as used in the handoff.

    For S69 this is ±(2,1), not the long complementary lobe.
    """
    a, b = deck_of_path(loop_fwd), deck_of_path(loop_rev)
    return a if abs(a[0]) + abs(a[1]) <= abs(b[0]) + abs(b[1]) else b


def preimage_record(blue, s, index):
    deg = akaho_degrees(s)
    loopA = branch_loop(blue, s, "A")
    loopB = branch_loop(blue, s, "B")
    dAB, dBA = deck_of_path(loopA), deck_of_path(loopB)
    return {
        "index": index,
        "pt": s["pt"],
        "kA": s["kA"],
        "tA": s["tA"],
        "kB": s["kB"],
        "tB": s["tB"],
        "dirA": s["dirA"],
        "dirB": s["dirB"],
        "phaseA": tangent_phase(s["dirA"]),
        "phaseB": tangent_phase(s["dirB"]),
        "degrees": deg,
        "deck_A_to_B": dAB,
        "deck_B_to_A": dBA,
        "deck_short": short_deck(loopA, loopB),
        "paramA": loc_param(s["kA"], s["tA"]),
        "paramB": loc_param(s["kB"], s["tB"]),
    }


def locate_orbit(orbs, target, tol=0.02):
    cands = sorted(
        ((i, pp, pre) for i, (pp, pre) in enumerate(orbs)),
        key=lambda t: _tdist(t[1], target),
    )
    i, pp, pre = cands[0]
    if _tdist(pp, target) > tol:
        return None
    return i, pp, pre


def build_q7(**kwargs):
    kw = dict(DEFAULT)
    kw.update(kwargs)
    red, blue, xinfo = build_pretzel(3, kw["blue_eps"], kw["red_eps"], kw["red_phi"])
    return red, blue, xinfo


def translate_to_kwz_chart(poly, special=KWZ_SPECIAL):
    """Translate a T^2 polyline so ``special`` becomes the chart point (0,0)."""
    return [((p[0] - special[0]) % TAU, (p[1] - special[1]) % TAU)
            for p in poly]


def translate_target_to_kwz(point, special=KWZ_SPECIAL):
    """Translate one pillowcase point into the admissible KWZ chart."""
    return P_point(((point[0] - special[0]) % TAU,
                    (point[1] - special[1]) % TAU))


def orbit_records(blue):
    orbs = orbit_group(blue)
    recs = []
    for i, (pp, pre) in enumerate(orbs):
        recs.append({
            "name": f"S{i}",
            "index": i,
            "P": tuple(pp),
            "n_pre": len(pre),
            "preimages": [preimage_record(blue, s, j) for j, s in enumerate(pre)],
        })
    return recs, orbs


def named_records(blue):
    recs, orbs = orbit_records(blue)
    out = {}
    for name, tgt in TARGETS.items():
        hit = locate_orbit(orbs, tgt)
        if hit is None:
            out[name] = None
            continue
        i, pp, pre = hit
        rec = recs[i]
        rec = dict(rec)
        rec["handoff_name"] = name
        rec["target"] = tgt
        rec["dist"] = _tdist(pp, tgt)
        out[name] = rec
    return out, recs, orbs


def find_crossing_by_edges(orbs, pair, tol=2):
    a, b = pair
    for i, (pp, pre) in enumerate(orbs):
        for s in pre:
            ka, kb = s["kA"], s["kB"]
            if abs(ka - a) <= tol and abs(kb - b) <= tol:
                return i, s
            if abs(ka - b) <= tol and abs(kb - a) <= tol:
                return i, s
    return None


def arc_passes_param(poly, k0, t0, k1, t1, forward, param, window=0.6):
    """Does the blue arc from (k0,t0) to (k1,t1) contain `param` as an interior
    point of the curve parameter k+t, without requiring a branch change?"""
    n = len(poly) - 1
    start = loc_param(k0, t0)
    end = loc_param(k1, t1)
    if forward:
        if end >= start:
            return start + 1e-9 < param < end - 1e-9
        return param > start + 1e-9 or param < end - 1e-9
    # backward
    if start >= end:
        return end + 1e-9 < param < start - 1e-9
    return param < start - 1e-9 or param > end + 1e-9


def accepted_self_bigon(blue, s0, s1, sel0, sel1, fb0, fb1):
    """One typed self-bigon: branch choices sel* in {0,1} and orientations fb*."""
    def halves(s, sel):
        if sel == 0:
            return (s["kA"], s["tA"]), (s["kB"], s["tB"])
        return (s["kB"], s["tB"]), (s["kA"], s["tA"])

    arr0, dep0 = halves(s0, sel0)
    arr1, dep1 = halves(s1, sel1)
    # same convention as maurer_cartan.self_polygon for k=2:
    # legs = (dep[0] -> arr[1], dep[1] -> arr[0])
    legs = [(dep0, arr1), (dep1, arr0)]
    fbs = (fb0, fb1)
    n = len(blue) - 1
    arcs = []
    for (p0, p1), fb in zip(legs, fbs):
        span = (p1[0] - p0[0]) % n if fb else (p0[0] - p1[0]) % n
        if span > 400:
            return None
        arcs.append(arc_between(blue, p0[0], p0[1], p1[0], p1[1], fb))
    loop, corners = _assemble_loop(arcs)
    if not bounds_disk(loop, corners):
        return None
    return {
        "legs": legs,
        "fbs": fbs,
        "sel": (sel0, sel1),
        "arcs": arcs,
    }


def find_handoff_witness(blue, orbs, s69_pre1):
    """Locate the straight-through self-bigon of handoff sec. 4."""
    hit0 = find_crossing_by_edges(orbs, WITNESS_EDGES[0])
    hit1 = find_crossing_by_edges(orbs, WITNESS_EDGES[1])
    if hit0 is None or hit1 is None:
        return None
    _, s0 = hit0
    _, s1 = hit1
    through = loc_param(s69_pre1["kB"], s69_pre1["tB"])
    # Handoff: branch states (1,0), orientations (forward, backward).
    typed = accepted_self_bigon(blue, s0, s1, 1, 0, True, False)
    if typed is None:
        # search all typings if the recorded one drifted
        for sel0 in (0, 1):
            for sel1 in (0, 1):
                for fb0 in (True, False):
                    for fb1 in (True, False):
                        typed = accepted_self_bigon(blue, s0, s1, sel0, sel1, fb0, fb1)
                        if typed is not None:
                            break
                    if typed is not None:
                        break
                if typed is not None:
                    break
            if typed is not None:
                break
    if typed is None:
        return dict(s0=s0, s1=s1, through=through, typed=None, hits=False)
    hits = []
    for i, ((p0, p1), fb) in enumerate(zip(typed["legs"], typed["fbs"])):
        if arc_passes_param(blue, p0[0], p0[1], p1[0], p1[1], fb, through):
            hits.append(i)
        # also test the recorded 941.640214
        if arc_passes_param(blue, p0[0], p0[1], p1[0], p1[1], fb, WITNESS_THROUGH_PARAM):
            if i not in hits:
                hits.append(i)
    return dict(s0=s0, s1=s1, through=through, typed=typed, hits=hits,
                idx0=hit0[0], idx1=hit1[0])


# ---------------------------------------------------------------------------
# Pillowcase-to-disk chart at the (0,0) corner (handoff sec. 7)
# ---------------------------------------------------------------------------
def r3_of(gamma, theta):
    return (math.cos(gamma), math.cos(theta), math.sin(gamma) * math.sin(theta))


def pillow_chart(gamma, theta, eps=1e-12):
    """Stereographic chart from the R^3 pillowcase model, centred at (0,0)."""
    v = r3_of(gamma, theta)
    nrm = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) or 1.0
    n = (v[0] / nrm, v[1] / nrm, v[2] / nrm)
    s2 = math.sqrt(2.0)
    n0 = (1.0 / s2, 1.0 / s2, 0.0)
    e1 = (1.0 / s2, -1.0 / s2, 0.0)
    e2 = (0.0, 0.0, 1.0)
    den = 1.0 - (n[0] * n0[0] + n[1] * n0[1] + n[2] * n0[2])
    if abs(den) < eps:
        return None
    x = (n[0] * e1[0] + n[1] * e1[1] + n[2] * e1[2]) / den
    y = (n[0] * e2[0] + n[1] * e2[1] + n[2] * e2[2]) / den
    return (x, y)


def chart_corner_images():
    """The four pillowcase corners in this chart.  (0,0) is the projection
    centre (sent to infinity).  The other three land on the x-axis."""
    return {
        (0.0, 0.0): pillow_chart(0.0, 0.0),
        (PI, 0.0): pillow_chart(PI, 0.0),
        (0.0, PI): pillow_chart(0.0, PI),
        (PI, PI): pillow_chart(PI, PI),
    }


def unwrap_seams(poly):
    """Lift a closed T^2 polyline across the seams gamma in {0, pi}.

    Returns a list of unwrapped R^2 points, so that a subsequent chart call
    can be applied after folding to a fundamental domain.  This is the
    bookkeeping step the handoff lists before the disk chart.
    """
    return unwrap_path(poly)


# ---------------------------------------------------------------------------
# KWZ skeleton encoding (handoff sec. 7, KWZ Def. 5.1--5.17)
# ---------------------------------------------------------------------------
# Disk chart: punctures at infinity (pillowcase (0,0)), (-1,0), (0,0), (1,0).
# Parameterizing arcs x = ±1/2 split the plane into three faces:
#   L : x < -1/2   (puncture (-1,0), D-loops)
#   M : |x| < 1/2  (puncture (0,0),  S-arrows)
#   R : x >  1/2   (puncture (1,0),  D-loops)
# Algebra words are pairs (face, winding) with multiplication
#   (f, m) * (g, n) = (f, m+n) if f == g else 0
# which is the handoff relation DS = SD = 0 together with D-powers
# (and S-powers) composing additively inside one face.

SKELETON = (-0.5, 0.5)
PUNCTURE = {"L": (-1.0, 0.0), "M": (0.0, 0.0), "R": (1.0, 0.0)}
CHART_JUMP = 0.75


def face_of(x, y, tol=1e-9):
    if x < -0.5 - tol:
        return "L"
    if x > 0.5 + tol:
        return "R"
    if abs(x) < 0.5 - tol:
        return "M"
    return None  # on a skeleton arc


def charted_polyline(blue, n_per_edge=6):
    """Dense chart image of the pillowcase curve, split at chart jumps.

    Each sample is (param, xy, Ppoint).  Returns a list of continuous
    segments (each a list of samples with xy != None).
    """
    n = len(blue) - 1
    raw = []
    for k in range(n):
        a = blue[k]
        b = blue[(k + 1) % n]
        db = ((b[0] - a[0] + PI) % TAU - PI, (b[1] - a[1] + PI) % TAU - PI)
        for s in range(n_per_edge):
            t = s / n_per_edge
            pt = ((a[0] + t * db[0]) % TAU, (a[1] + t * db[1]) % TAU)
            pp = P_point(pt)
            xy = pillow_chart(pp[0], pp[1])
            raw.append((k + t, xy, pp))
    segs, cur = [], []
    for i, samp in enumerate(raw):
        xy = samp[1]
        if xy is None:
            if cur:
                segs.append(cur)
                cur = []
            continue
        if cur:
            prev = cur[-1][1]
            if math.hypot(xy[0] - prev[0], xy[1] - prev[1]) > CHART_JUMP:
                segs.append(cur)
                cur = []
        cur.append(samp)
    if cur:
        segs.append(cur)
    # If the first and last segments meet continuously, glue them (closed curve).
    if len(segs) >= 2:
        a, b = segs[0][0][1], segs[-1][-1][1]
        if math.hypot(a[0] - b[0], a[1] - b[1]) < CHART_JUMP:
            segs[0] = segs[-1] + segs[0]
            segs.pop()
    return segs


def _cross_x(p, q, xc):
    """If the segment p->q in the chart crosses the vertical line x=xc,
    return (param, y, direction).  direction is +1 if x increases through xc."""
    (tp, (xp, yp), _), (tq, (xq, yq), _) = p, q
    if (xp - xc) * (xq - xc) >= 0:
        return None
    lam = (xc - xp) / (xq - xp)
    y = yp + lam * (yq - yp)
    return (tp + lam * (tq - tp), y, 1 if xq > xp else -1)


def skeleton_hits(segs):
    """Legacy chart-segment finder, retained only for comparison.

    This misses intersections across ``CHART_JUMP`` breaks and must not be
    used to construct the type-D object.  Use :func:`torus_skeleton_hits`.
    """
    hits = []
    for seg in segs:
        for i in range(len(seg) - 1):
            for xc, arc in ((-0.5, "L"), (0.5, "R")):
                cr = _cross_x(seg[i], seg[i + 1], xc)
                if cr is None:
                    continue
                param, y, dirc = cr
                hits.append({
                    "arc": arc,
                    "x": xc,
                    "y": y,
                    "param": param,
                    "dir": dirc,  # +1: left->right through the arc
                    "pt": (xc, y),
                })
    hits.sort(key=lambda h: h["param"])
    # drop near-duplicates from densification
    out = []
    for h in hits:
        if out and abs(h["param"] - out[-1]["param"]) < 1e-4 and h["arc"] == out[-1]["arc"]:
            continue
        out.append(h)
    return out


def _edge_point(blue, k, t):
    """Point at parameter ``k+t`` using the short lift of a T^2 edge."""
    n = len(blue) - 1
    a, b = blue[k], blue[(k + 1) % n]
    dg = (b[0] - a[0] + PI) % TAU - PI
    dt = (b[1] - a[1] + PI) % TAU - PI
    return ((a[0] + t * dg) % TAU, (a[1] + t * dt) % TAU)


def _edge_chart(blue, k, t):
    pt = _edge_point(blue, k, t)
    return pillow_chart(*P_point(pt))


def _bisect_skeleton(blue, k, ta, tb, xc, max_abs=float("inf")):
    """Bisect one genuine finite crossing of ``x=xc`` on a T^2 edge.

    A sign change caused by running through the stereographic projection
    centre is not a crossing of the parameterizing arc.  Such a bracket hits
    the undefined chart point and is rejected.  ``max_abs`` is retained only
    as an optional diagnostic; the construction itself uses no finite cutoff.
    """
    xa = _edge_chart(blue, k, ta)
    xb = _edge_chart(blue, k, tb)
    if xa is None or xb is None:
        return None
    fa, fb = xa[0] - xc, xb[0] - xc
    if fa == 0.0:
        tm, xy = ta, xa
    elif fb == 0.0:
        tm, xy = tb, xb
    elif fa * fb > 0.0:
        return None
    else:
        for _ in range(55):
            tm = 0.5 * (ta + tb)
            xy = _edge_chart(blue, k, tm)
            if xy is None or max(abs(xy[0]), abs(xy[1])) > max_abs:
                return None
            fm = xy[0] - xc
            if abs(fm) < 1e-13:
                break
            if fa * fm <= 0.0:
                tb, xb, fb = tm, xy, fm
            else:
                ta, xa, fa = tm, xy, fm
        tm = 0.5 * (ta + tb)
        xy = _edge_chart(blue, k, tm)
    if xy is None or max(abs(xy[0]), abs(xy[1])) > max_abs:
        return None
    # A discontinuity at the projection centre can also produce a sign-changing
    # bracket whose bisection limit is a pole rather than a root.  Accept only
    # an actual finite solution of x=xc.
    if abs(xy[0] - xc) > 1e-8:
        return None
    # Direction in the chart, evaluated on a stable edge-local bracket.  Do
    # not tie this probe to the final bisection bracket: after 55 iterations
    # that bracket is at roundoff scale and can reverse a transverse sign.
    eps = min(1e-5, 0.25 * tm, 0.25 * (1.0 - tm))
    if eps <= 1e-12:
        eps = 1e-8
    lo = _edge_chart(blue, k, max(0.0, tm - eps))
    hi = _edge_chart(blue, k, min(1.0, tm + eps))
    if lo is None or hi is None:
        return None
    return tm, xy[1], 1 if hi[0] > lo[0] else -1


def torus_skeleton_hits(blue, n_samp=40, max_abs=float("inf")):
    """Find x=±1/2 intersections directly on the original T^2 edges.

    Subdivision is used only to bracket roots; each accepted intersection is
    then bisected on the exact PL edge.  This makes crossings independent of
    the chart-space jump splitting that corrupted the earlier 25-dot count.
    The default has no affine-chart cutoff: finite roots high in the chart are
    genuine skeleton intersections, even near the deleted projection corner.
    """
    hits = []
    n = len(blue) - 1
    for k in range(n):
        vals = []
        for s in range(n_samp + 1):
            t = s / n_samp
            vals.append((t, _edge_chart(blue, k, t)))
        for s in range(n_samp):
            ta, a = vals[s]
            tb, b = vals[s + 1]
            if a is None or b is None:
                continue
            for xc, arc in ((-0.5, "L"), (0.5, "R")):
                if (a[0] - xc) * (b[0] - xc) > 0.0:
                    continue
                root = _bisect_skeleton(blue, k, ta, tb, xc, max_abs=max_abs)
                if root is None:
                    continue
                t, y, direction = root
                h = {"arc": arc, "x": xc, "y": y, "param": k + t,
                     "dir": direction, "pt": (xc, y)}
                if hits and hits[-1]["arc"] == arc and abs(hits[-1]["param"] - h["param"]) < 1e-8:
                    continue
                hits.append(h)
    hits.sort(key=lambda h: h["param"])
    return hits


def angular_change(path_xy, puncture):
    """Unwrapped angular change of a chart path about ``puncture``."""
    if len(path_xy) < 2:
        return 0.0
    ang = 0.0
    px, py = puncture
    prev = math.atan2(path_xy[0][1] - py, path_xy[0][0] - px)
    for x, y in path_xy[1:]:
        a = math.atan2(y - py, x - px)
        dang = (a - prev + PI) % TAU - PI
        ang += dang
        prev = a
    return ang


def winding_about(path_xy, puncture):
    """Net full winding of a chart path about a puncture."""
    return int(round(angular_change(path_xy, puncture) / TAU))


def _vertical_segment(x, y0, y1, n=16):
    return [(x, y0 + (y1 - y0) * k / n) for k in range(n + 1)]


def standardized_face_turn(path_xy, start_arc, end_arc, face):
    """Relative angular change after standardizing skeleton endpoints.

    Each vertical skeleton side is based at ``y=0``.  Adding the two side
    segments turns an open face join into a path between the standard quiver
    vertices, whose angular change is an integral multiple of ``2*pi/n_f``.
    """
    if len(path_xy) < 2:
        return 0.0
    x0 = -0.5 if start_arc == "L" else 0.5
    x1 = -0.5 if end_arc == "L" else 0.5
    y0, y1 = path_xy[0][1], path_xy[-1][1]
    closed = (_vertical_segment(x0, 0.0, y0)[:-1] + path_xy
              + _vertical_segment(x1, y1, 0.0)[1:])
    return angular_change(closed, PUNCTURE[face])


def face_path_word(path_xy, start_arc, end_arc, face):
    """Cyclic-quiver word ``(face, length)`` for an oriented face join.

    Outer faces have one side, hence one arrow per full turn.  The middle
    face has two sides, hence one arrow per half turn.  Length zero is an
    identity component and negative length is clockwise.
    """
    turn = standardized_face_turn(path_xy, start_arc, end_arc, face)
    unit = TAU if face in ("L", "R") else PI
    length = int(round(turn / unit))
    return (face, length)


def segment_word(seg_xy, start_arc, end_arc, start_dir=None):
    """Algebra word for a chart path between two skeleton hits.

    Returns ``(face, cyclic-quiver length)``.  Positive length is
    anticlockwise, negative length clockwise, and zero is an identity
    component removed by geometric bigon cancellation.
    """
    if not seg_xy:
        return None
    # Consecutive hits on distinct arcs necessarily bound a middle-face
    # segment.  For two hits on the same arc, the direction at the first hit
    # says on which side the segment departs.  This topological rule remains
    # valid when an outer-face path crosses the chart cut at infinity.
    if start_arc != end_arc:
        face = "M"
    elif start_arc == "L":
        face = "L" if start_dir == -1 else "M"
    else:
        face = "R" if start_dir == 1 else "M"
    return face_path_word(seg_xy, start_arc, end_arc, face)


def torus_chart_path_between(blue, t0, t1, subdivisions=6):
    """Chart the oriented T^2 path from curve parameter ``t0`` to ``t1``.

    Unlike ``chart_path_between``, this does not discard the two sides of a
    stereographic chart jump.  Angles used by ``winding_about`` are unwrapped
    across that jump, so a finite path around the puncture at infinity remains
    available for the outer-face word.
    """
    n = len(blue) - 1
    if t1 <= t0:
        t1 += n
    cuts = [t0]
    first_edge = int(math.floor(t0)) + 1
    last_edge = int(math.floor(t1))
    cuts.extend(float(k) for k in range(first_edge, last_edge + 1))
    cuts.append(t1)
    pts = []
    for a, b in zip(cuts, cuts[1:]):
        if b - a < 1e-12:
            continue
        for s in range(subdivisions):
            p = a + (b - a) * s / subdivisions
            k = int(math.floor(p)) % n
            u = p - math.floor(p)
            xy = _edge_chart(blue, k, u)
            if xy is not None:
                pts.append(xy)
    p = t1
    k = int(math.floor(p)) % n
    u = p - math.floor(p)
    xy = _edge_chart(blue, k, u)
    if xy is not None:
        pts.append(xy)
    return pts


def chart_path_between(segs, t0, t1):
    """Chart samples with param in (t0, t1) along the (possibly wrapped) curve."""
    pts = []
    for seg in segs:
        for samp in seg:
            t, xy, _ = samp
            if xy is None:
                continue
            if t0 <= t1:
                if t0 < t < t1:
                    pts.append(xy)
            else:
                if t > t0 or t < t1:
                    pts.append(xy)
    return pts


def unique_hits(hits, ytol=0.001):
    """Geometric skeleton dots, forgetting the double cover of P by T^2."""
    out = []
    for h in hits:
        found = False
        for u in out:
            if h["arc"] == u["arc"] and abs(h["y"] - u["y"]) < ytol:
                u["params"].append(h["param"])
                found = True
                break
        if not found:
            rec = dict(h)
            rec["params"] = [h["param"]]
            out.append(rec)
    return out


def encode_type_d(blue, n_samp=40):
    """Build the KWZ precurve of the q=7 blue curve.

    Returns dict with hits, oriented words along the curve, and the sparse
    maps delta_fwd / delta_rev : i -> (word, j).
    """
    segs = charted_polyline(blue)
    hits = torus_skeleton_hits(blue, n_samp=n_samp)
    words = []
    for i in range(len(hits)):
        j = (i + 1) % len(hits)
        path = torus_chart_path_between(blue, hits[i]["param"], hits[j]["param"])
        w = segment_word(path, hits[i]["arc"], hits[j]["arc"], hits[i]["dir"])
        words.append(w)
    return {
        "segs": segs,
        "hits": hits,
        "unique": unique_hits(hits),
        "words": words,
    }


def mul_word(a, b):
    """Cyclic-quiver multiplication: concatenate within a common face."""
    if a is None or b is None:
        return None
    if a[0] == "1":
        return b
    if b[0] == "1":
        return a
    if a[0] != b[0]:
        return None
    return (a[0], a[1] + b[1])


def delta_squared_words(words):
    """The two-step compositions that would appear in δ² along the oriented curve."""
    n = len(words)
    out = []
    for i in range(n):
        prod = mul_word(words[i], words[(i + 1) % n])
        if prod is not None and prod[0] != "1":
            out.append((i, prod))
    return out


def _match_unique(h, uniq, ytol=0.001):
    for i, u in enumerate(uniq):
        if h["arc"] == u["arc"] and abs(h["y"] - u["y"]) < ytol:
            return i
    return None


def single_copy_precurve(data, ytol=0.001):
    """Legacy first-return attempt, retained as a negative diagnostic.

    The T^2 walk covers P twice.  Walking the skeleton hits in parameter
    order and recording first returns to the starting geometric dot gives
    a single cycle on the unique dots (or a diagnostic if it does not).
    """
    hits, uniq = data["hits"], data["unique"]
    seq = []
    for h in hits:
        i = _match_unique(h, uniq, ytol)
        if i is None:
            continue
        if not seq or seq[-1] != i:
            seq.append(i)
    if len(seq) >= 2 and seq[0] in seq[1:]:
        cut = seq[1:].index(seq[0]) + 1
        cycle = seq[:cut]
    else:
        cycle = seq
    # Drop a closing repeat if the walk recorded the start twice.
    if len(cycle) >= 2 and cycle[0] == cycle[-1]:
        cycle = cycle[:-1]
    # If one extra vertex is a duplicate interior, keep the first occurrence
    # of each index in cycle order (a single traversal cannot repeat a dot).
    seen, dedup = set(), []
    for i in cycle:
        if i in seen:
            continue
        seen.add(i)
        dedup.append(i)
    extras = [i for i in cycle if cycle.count(i) > 1]
    cycle = dedup
    # words along this cycle: take the T^2 word at the first time each
    # consecutive pair appears in the hit list
    pair_word = {}
    for a, b, w in zip(range(len(hits)), range(1, len(hits) + 1), data["words"]):
        ia = _match_unique(hits[a], uniq, ytol)
        ib = _match_unique(hits[b % len(hits)], uniq, ytol)
        if ia is None or ib is None or ia == ib:
            continue
        pair_word.setdefault((ia, ib), w)
    words = []
    missing = []
    for k in range(len(cycle)):
        ia, ib = cycle[k], cycle[(k + 1) % len(cycle)]
        w = pair_word.get((ia, ib))
        words.append(w)
        if w is None:
            missing.append((ia, ib))
    # δ as a list of (src, word, tgt) along the oriented cycle
    _ = extras  # used in the report via the returned dict
    delta = []
    for k, w in enumerate(words):
        if w is None:
            continue
        delta.append((cycle[k], w, cycle[(k + 1) % len(cycle)]))
    return {
        "uniq": uniq,
        "cycle": cycle,
        "words": words,
        "delta": delta,
        "missing": missing,
        "seq": seq,
    }


def geometric_precurve(data, ytol=0.001):
    """Construct and fully cancel the precurve on the pillowcase.

    The closed T^2 parameter walk traverses every geometric face join twice,
    in opposite directions.  After projection its non-loop join graph is an
    interval; the two loops at its degree-one vertices run into deleted
    pillowcase corners and are one-sided joins.  A zero-length two-sided join
    is a removable bigon with the skeleton.  Cancelling such a bigon removes
    its two dots and concatenates the adjacent joins in their common face.
    """
    hits, uniq, words = data["hits"], data["unique"], data["words"]
    seq = [_match_unique(h, uniq, ytol) for h in hits]
    groups = {}
    for k in range(len(hits)):
        i, j = seq[k], seq[(k + 1) % len(hits)]
        key = (min(i, j), max(i, j))
        groups.setdefault(key, []).append((i, words[k], j, k))

    loops = {key[0]: vals for key, vals in groups.items() if key[0] == key[1]}
    joins = {key: vals for key, vals in groups.items() if key[0] != key[1]}
    adjacency = {i: [] for i in range(len(uniq))}
    for i, j in joins:
        adjacency[i].append(j)
        adjacency[j].append(i)
    endpoints = sorted(i for i, nbrs in adjacency.items() if len(nbrs) == 1)

    path = []
    if len(endpoints) == 2:
        prev, cur = None, endpoints[0]
        while True:
            path.append(cur)
            nxt = [j for j in adjacency[cur] if j != prev]
            if not nxt:
                break
            if len(nxt) != 1:
                break
            prev, cur = cur, nxt[0]

    signed_words = []
    paired_ok = True
    for i, j in zip(path, path[1:]):
        vals = joins.get((min(i, j), max(i, j)), [])
        fwd = [w for a, w, b, _ in vals if a == i and b == j]
        rev = [w for a, w, b, _ in vals if a == j and b == i]
        if len(fwd) != 1 or len(rev) != 1:
            paired_ok = False
            signed_words.append(None)
            continue
        a, b = fwd[0], rev[0]
        if a[0] != b[0] or a[1] != -b[1]:
            paired_ok = False
        signed_words.append(a)

    vertices = list(path)
    edges = list(signed_words)
    cancellations = []
    while edges and any(w is not None and w[1] == 0 for w in edges):
        candidates = [k for k, w in enumerate(edges)
                      if w is not None and w[1] == 0
                      and 0 < k and k + 1 < len(vertices) - 1]
        if not candidates:
            break
        k = candidates[0]
        left, zero, right = edges[k - 1], edges[k], edges[k + 1]
        if left is None or right is None or left[0] != right[0]:
            break
        cancellations.append({
            "dots": (vertices[k], vertices[k + 1]),
            "zero_face": zero[0],
            "joined_face": left[0],
            "joined_length": left[1] + right[1],
        })
        new_word = (left[0], left[1] + right[1])
        vertices = vertices[:k] + vertices[k + 2:]
        edges = edges[:k - 1] + [new_word] + edges[k + 2:]

    delta = []
    for i, w, j in zip(vertices, edges, vertices[1:]):
        if w is None or w[1] == 0:
            continue
        if w[1] > 0:
            delta.append((i, w, j))
        else:
            delta.append((j, (w[0], -w[1]), i))
    delta = collect_f2(delta)
    return {
        "uniq": uniq,
        "seq": seq,
        "groups": groups,
        "loops": loops,
        "joins": joins,
        "adjacency": adjacency,
        "endpoints": endpoints,
        "path": path,
        "signed_words": signed_words,
        "paired_ok": paired_ok,
        "vertices": vertices,
        "words": edges,
        "cancellations": cancellations,
        "delta": delta,
        "residue": apply_mul(delta, delta),
    }


def cyclic_precurve(data):
    """Cancel zero joins in a single oriented cyclic precurve.

    Unlike :func:`geometric_precurve`, this applies when one T^2 component
    projects one-to-one to a pillowcase circle.  A zero join removes its two
    endpoint dots; the adjacent words concatenate in their common face.
    Original dot labels are retained so that two cyclic words can be compared
    directly.
    """
    n = len(data["unique"])
    if len(data["hits"]) != n or len(data["words"]) != n:
        raise AssertionError("cyclic precurve requires one hit per geometric dot")
    successor = {
        i: (data["words"][i], (i + 1) % n)
        for i in range(n)
    }
    cancellations = []
    while True:
        zeros = [i for i, (word, _) in successor.items() if word[1] == 0]
        if not zeros:
            break
        src = zeros[0]
        _, tgt = successor[src]
        pred = next(i for i, (_, j) in successor.items() if j == src)
        left, _ = successor[pred]
        right, out = successor[tgt]
        if left[0] != right[0]:
            raise AssertionError("zero cancellation has incompatible adjacent faces")
        joined = (left[0], left[1] + right[1])
        cancellations.append({
            "dots": (src, tgt),
            "zero_face": successor[src][0][0],
            "joined_face": left[0],
            "joined_length": joined[1],
        })
        del successor[src]
        del successor[tgt]
        successor[pred] = (joined, out)

    start = min(successor)
    vertices, words = [], []
    cur = start
    while True:
        vertices.append(cur)
        word, nxt = successor[cur]
        words.append(word)
        cur = nxt
        if cur == start:
            break
        if len(vertices) > len(successor):
            raise AssertionError("cancelled cyclic join graph is not one cycle")
    if len(vertices) != len(successor):
        raise AssertionError("cancelled cyclic join graph is disconnected")

    delta = []
    for src, word, tgt in zip(vertices, words, vertices[1:] + vertices[:1]):
        if word[1] > 0:
            delta.append((src, word, tgt))
        elif word[1] < 0:
            delta.append((tgt, (word[0], -word[1]), src))
        else:
            raise AssertionError("zero word survived cyclic cancellation")
    delta = collect_f2(delta)
    return {
        "vertices": vertices,
        "words": words,
        "cancellations": cancellations,
        "delta": delta,
        "residue": apply_mul(delta, delta),
    }


def collect_f2(triples):
    acc = {}
    for i, w, j in triples:
        key = (i, w, j)
        acc[key] = acc.get(key, 0) ^ 1
    return [(i, w, j) for (i, w, j), v in acc.items() if v]


def apply_mul(left, right):
    """Compose two morphisms (lists of (i, word, j)) over F2."""
    out = []
    for i, a, k in left:
        for kk, b, j in right:
            if k != kk:
                continue
            p = mul_word(a, b)
            if p is not None:
                out.append((i, p, j))
    return collect_f2(out)


def nearest_dot(uniq, xy, arc=None, ytol=0.25):
    cands = [(i, u) for i, u in enumerate(uniq)
             if arc is None or u["arc"] == arc]
    if not cands:
        return None
    i, u = min(cands, key=lambda t: abs(t[1]["y"] - xy[1]))
    if abs(u["y"] - xy[1]) > ytol:
        return None
    return i, u


def resolve_to_morphism(blue, pre, face, xy, puncture, xc, uniq):
    """φ(x) as a list of (src, word, tgt) on the unique-dot index set."""
    res = resolve_crossing_in_face(blue, pre, face, xy, puncture, xc)
    terms = []
    unmatched = []
    for tag, wind, a, b in res.get("terms", []):
        if math.hypot(a[0] - b[0], a[1] - b[1]) < 0.05:
            continue  # identity / Mor^×
        ia = nearest_dot(uniq, a, arc="R" if face == "R" else "L")
        ib = nearest_dot(uniq, b, arc="R" if face == "R" else "L")
        if ia is None or ib is None:
            unmatched.append((tag, wind, a, b))
            continue
        # face-local word: D-power or S-power of this winding
        word = (face, wind)
        terms.append((ia[0], word, ib[0]))
    return {
        "res": res,
        "u": collect_f2(terms),
        "unmatched": unmatched,
    }


def mc_residue(delta, u):
    """D_N(u) + u^2 over F2, as a list of (i, word, j)."""
    du = collect_f2(apply_mul(delta, u) + apply_mul(u, delta))
    uu = apply_mul(u, u)
    return collect_f2(du + uu), du, uu


def word_bigrading(word):
    """KWZ (quantum, delta) degree of a positive cyclic-quiver word."""
    face, length = word
    qdeg = Fraction(-2 * length, 1 if face in ("L", "R") else 2)
    return qdeg, qdeg / 2


def generator_bigradings(vertices, delta):
    """Solve relative generator bigradings from the type-D differential."""
    adjacency = {i: [] for i in vertices}
    for i, word, j in delta:
        qword, dword = word_bigrading(word)
        # q(j)=q(i)-q(word), delta(j)=delta(i)-1-delta(word)
        adjacency[i].append((j, -qword, -1 - dword))
        adjacency[j].append((i, qword, 1 + dword))
    grades = {vertices[0]: (Fraction(0), Fraction(0))}
    stack = [vertices[0]]
    consistent = True
    while stack:
        i = stack.pop()
        qi, di = grades[i]
        for j, dq, dd in adjacency[i]:
            value = (qi + dq, di + dd)
            if j in grades:
                if grades[j] != value:
                    consistent = False
            else:
                grades[j] = value
                stack.append(j)
    return grades, consistent


def morphism_term_bigrading(term, grades):
    i, word, j = term
    qw, dw = word_bigrading(word)
    qi, di = grades[i]
    qj, dj = grades[j]
    return qw + qj - qi, dw + dj - di


def end_terms_in_delta_degree(data, pre, grades, delta_degree):
    """All valid End terms of one fixed relative delta degree.

    Fixing delta degree determines the positive path length from the endpoint
    gradings, so this basis is finite even though the face algebras are not.
    """
    delta_degree = Fraction(delta_degree)
    terms = []
    for i in pre["vertices"]:
        arc_i = data["unique"][i]["arc"]
        side_i = 0 if arc_i == "L" else 1
        qi, di = grades[i]
        for j in pre["vertices"]:
            arc_j = data["unique"][j]["arc"]
            side_j = 0 if arc_j == "L" else 1
            qj, dj = grades[j]
            if arc_i == arc_j and dj - di == delta_degree:
                terms.append((i, ("1", 0), j))
            for face in ("L", "M", "R"):
                if face == "L" and (arc_i != "L" or arc_j != "L"):
                    continue
                if face == "R" and (arc_i != "R" or arc_j != "R"):
                    continue
                scale = 1 if face in ("L", "R") else 2
                length = scale * (dj - di - delta_degree)
                if length.denominator != 1 or length <= 0:
                    continue
                length = int(length)
                if face == "M" and (side_i + length - side_j) % 2:
                    continue
                terms.append((i, (face, length), j))
    return terms


def end_differential(delta, terms):
    return collect_f2(apply_mul(delta, terms) + apply_mul(terms, delta))


def sparse_column_rank(columns):
    """Rank over F2 of columns given as iterables of row indices."""
    pivots = {}
    rank = 0
    for rows in columns:
        value = 0
        for row in rows:
            value ^= 1 << row
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                rank += 1
                break
    return rank


def end_h_representatives():
    """A homogeneous basis for H^(q,-1) End(N), keyed by quantum degree."""
    return {
        Fraction(-16): [[(19, ("L", 1), 37)]],
        Fraction(-14): [[(7, ("M", 2), 32), (4, ("M", 2), 33)]],
        Fraction(-12): [[(7, ("M", 2), 27), (4, ("M", 2), 24)]],
        Fraction(-10): [[(1, ("M", 1), 32), (0, ("M", 2), 33)]],
        Fraction(-8): [
            [(1, ("M", 1), 33)],
            [(1, ("M", 1), 27), (0, ("M", 2), 24)],
        ],
        Fraction(-6): [[(1, ("M", 1), 24)]],
        Fraction(-2): [
            [(19, ("L", 1), 19)],
            [(7, ("M", 2), 7), (4, ("M", 2), 4)],
            [(22, ("M", 1), 24)],
            [(32, ("M", 2), 32), (33, ("M", 2), 33)],
            [(37, ("L", 1), 37)],
        ],
        Fraction(0): [[(7, ("M", 2), 12), (4, ("M", 2), 15)]],
        Fraction(2): [
            [(7, ("M", 2), 17), (4, ("M", 2), 16)],
            [(1, ("M", 1), 7), (0, ("M", 2), 4)],
        ],
        Fraction(4): [[(1, ("M", 1), 12), (0, ("M", 2), 15)]],
        Fraction(6): [
            [(1, ("M", 1), 15)],
            [(1, ("M", 1), 17), (0, ("M", 2), 16)],
        ],
        Fraction(8): [[(1, ("M", 1), 16)]],
        Fraction(10): [
            [(22, ("M", 1), 15)],
            [(32, ("M", 2), 7), (33, ("M", 2), 4)],
        ],
        Fraction(12): [[(37, ("L", 1), 19)]],
    }


def report_encoding(blue, s69_target=None):
    print("== KWZ type-D encoding of the q=7 blue curve ==")
    data = encode_type_d(blue)
    hits, words, uniq = data["hits"], data["words"], data["unique"]
    # Root locations, rather than subdivision samples, are the output.  The
    # three independent bracket resolutions must therefore give the same
    # ordered intersection record.
    stable = []
    for ns in (20, 40, 80):
        hs = torus_skeleton_hits(blue, n_samp=ns)
        us = unique_hits(hs)
        stable.append((ns, hs, us))
    sigs = [tuple((h["arc"], round(h["param"], 8), round(h["y"], 7))
                  for h in hs) for _, hs, _ in stable]
    print("  subdivision stability: " + ", ".join(
        f"n={ns}: {len(hs)} hits/{len(us)} dots" for ns, hs, us in stable))
    check(sigs[0] == sigs[1] == sigs[2],
          "T^2 skeleton hits are stable for n_samp in {20,40,80}")
    multiplicities = [len(u["params"]) for u in uniq]
    check(multiplicities and all(m == 2 for m in multiplicities),
          "every geometric skeleton dot has multiplicity exactly two")
    word_sigs = []
    for ns in (3, 6, 12, 24):
        ws = []
        for k in range(len(hits)):
            j = (k + 1) % len(hits)
            path = torus_chart_path_between(
                blue, hits[k]["param"], hits[j]["param"], subdivisions=ns)
            ws.append(segment_word(
                path, hits[k]["arc"], hits[j]["arc"], hits[k]["dir"]))
        word_sigs.append(ws)
    check(all(ws == words for ws in word_sigs),
          "all cyclic-quiver words are stable at path subdivisions 3,6,12,24")
    print(f"  charted segments: {len(data['segs'])}")
    print(f"  skeleton intersections (T^2 walk): {len(hits)}")
    nL = sum(1 for h in hits if h["arc"] == "L")
    nR = sum(1 for h in hits if h["arc"] == "R")
    print(f"    on x=-1/2: {nL}    on x=+1/2: {nR}")
    uL = sum(1 for h in uniq if h["arc"] == "L")
    uR = sum(1 for h in uniq if h["arc"] == "R")
    print(f"  geometric dots on P (deduped): {len(uniq)}  "
          f"(L {uL}, R {uR})")
    dbl = [u for u in uniq if len(u["params"]) >= 2]
    print(f"    dots seen twice on the T^2 walk: {len(dbl)}")
    from collections import Counter
    wc = Counter(words)
    print(f"  oriented face words: {dict(wc)}")
    check(len(hits) > 0, "the charted curve meets the skeleton")
    rec = named_s69_chart(blue, target=s69_target)
    if rec is not None:
        print(f"  S69 in the disk chart: xy={tuple(round(v, 5) for v in rec['xy'])}  "
              f"face={rec['face']}")
        check(rec["face"] is not None, "S69 is not on the skeleton")
    pre = geometric_precurve(data)
    print(f"  projected join graph: {len(pre['joins'])} two-sided joins, "
          f"{len(pre['loops'])} one-sided endpoint joins")
    print(f"    endpoints={pre['endpoints']}  path length={len(pre['path'])}")
    check(len(pre["joins"]) == 36 and len(pre["loops"]) == 2,
          "the projected graph has 36 two-sided and two endpoint joins")
    check(pre["paired_ok"],
          "every two-sided join has two opposite T^2 orientations")
    check(len(pre["path"]) == len(uniq)
          and set(pre["path"]) == set(range(len(uniq))),
          "the non-loop geometric join graph is a 37-dot interval")
    check(set(pre["loops"]) == set(pre["endpoints"]),
          "the one-sided joins occur exactly at the interval endpoints")
    print(f"  zero-length bigon cancellations: {pre['cancellations']}")
    check(len(pre["cancellations"]) == 3
          and len(pre["vertices"]) == 31
          and all(w is not None and w[1] != 0 for w in pre["words"]),
          "three bigon cancellations give a fully cancelled 31-dot precurve")
    print(f"  fully cancelled delta: |V|={len(pre['vertices'])}, "
          f"|delta|={len(pre['delta'])}, |delta^2|={len(pre['residue'])}")
    check(len(pre["delta"]) == 30 and not pre["residue"],
          "the 31-dot type-D differential satisfies delta^2=0")
    data["precurve"] = pre
    return data


def named_s69_chart(blue, target=None):
    recs, orbs = orbit_records(blue)
    hit = locate_orbit(orbs, TARGETS["S69"] if target is None else target)
    if hit is None:
        return None
    _, pp, _ = hit
    xy = pillow_chart(pp[0], pp[1])
    if xy is None:
        return {"xy": None, "face": None}
    return {"xy": xy, "face": face_of(xy[0], xy[1]), "P": pp}


def face_arcs(segs, face, xc):
    """Maximal charted subpaths that stay in `face` and meet x=xc at both ends.

    Each returned arc is a list of (x,y) including the two skeleton endpoints.
    """
    out = []
    for seg in segs:
        i = 0
        while i < len(seg) - 1:
            xy = seg[i][1]
            if face_of(xy[0], xy[1]) != face:
                i += 1
                continue
            j = i
            while j < len(seg) - 1 and face_of(seg[j + 1][1][0], seg[j + 1][1][1]) == face:
                j += 1
            pts = [seg[k][1] for k in range(i, j + 1)]
            if len(pts) >= 2:
                # snap ends to the skeleton if they are close
                if abs(pts[0][0] - xc) < 0.15 and abs(pts[-1][0] - xc) < 0.15:
                    out.append(pts)
            i = j + 1
    return out


def closest_on_arc(arc, xy):
    best, bi = 1e9, 0
    for i, p in enumerate(arc):
        d = math.hypot(p[0] - xy[0], p[1] - xy[1])
        if d < best:
            best, bi = d, i
    return bi, best


def right_turn_reconnect(arc1, i1, arc2, i2):
    """The two right-turn smoothings of two arcs crossing at indices i1, i2.

    A right turn takes the incoming tangent of arc1 and rotates clockwise onto
    arc2.  We return both possible pairings as chart paths; KWZ keeps the
    anticlockwise ones about the face puncture.
    """
    # pairing 1: head of arc1 + tail of arc2, and vice versa
    a = list(reversed(arc1[:i1 + 1])) + arc2[i2:]
    b = list(reversed(arc2[:i2 + 1])) + arc1[i1:]
    # pairing 2
    c = arc1[:i1 + 1] + arc2[i2:]
    d = arc2[:i2 + 1] + arc1[i1:]
    return [("pair1", a), ("pair1", b), ("pair2", c), ("pair2", d)]


def _pt_on_blue(blue, k, t):
    n = len(blue) - 1
    a = blue[k]
    b = blue[(k + 1) % n]
    db = ((b[0] - a[0] + PI) % TAU - PI, (b[1] - a[1] + PI) % TAU - PI)
    return ((a[0] + t * db[0]) % TAU, (a[1] + t * db[1]) % TAU)


def walk_branch_to_skeleton(blue, k, t, forward, face, xc, max_steps=4000):
    """Walk the T^2 curve from (k, t) until the chart hits x=xc or leaves `face`."""
    n = len(blue) - 1
    path = []
    step = 1 if forward else -1
    cur_k, cur_t = k, t
    for _ in range(max_steps):
        pt = _pt_on_blue(blue, cur_k, cur_t)
        xy = pillow_chart(*P_point(pt))
        if xy is None:
            return None
        path.append(xy)
        fc = face_of(xy[0], xy[1])
        if len(path) > 3 and abs(xy[0] - xc) < 0.08:
            return path
        if len(path) > 8 and fc is not None and fc != face:
            return None
        # advance ~0.15 of an edge
        if forward:
            cur_t += 0.2
            if cur_t >= 1.0:
                cur_t -= 1.0
                cur_k = (cur_k + 1) % n
        else:
            cur_t -= 0.2
            if cur_t <= 0.0:
                cur_t += 1.0
                cur_k = (cur_k - 1) % n
    return None


def resolve_crossing_in_face(blue, s_preimages, face, xy, puncture, xc):
    """KWZ Def. 5.17 resolution of a geometric crossing that sits in `face`.

    From each half-edge at the two T^2 preimages, walk in both directions
    until the chart hits the skeleton, staying in `face`.  Pair those rays
    by the right-turn rule and keep anticlockwise face words.
    """
    rays = []
    for s in s_preimages:
        for br, kk, tt in (("A", s["kA"], s["tA"]), ("B", s["kB"], s["tB"])):
            for fwd in (True, False):
                path = walk_branch_to_skeleton(blue, kk, tt, fwd, face, xc)
                if path is not None and len(path) >= 2:
                    rays.append({
                        "branch": br,
                        "fwd": fwd,
                        "path": path,
                        "end": path[-1],
                        "pt": s["pt"],
                    })
    # four geometric half-edges live on P; the T^2 double cover may duplicate them
    print(f"  rays that stay in face {face} and reach x={xc}: {len(rays)}")
    for r in rays:
        print(f"    {r['branch']} {'fwd' if r['fwd'] else 'bwd'}  "
              f"len={len(r['path'])}  end={tuple(round(v, 3) for v in r['end'])}")
    if len(rays) < 4:
        return {"carriers": len(rays) // 2, "terms": [], "rays": rays}
    # Right-turn pairings: incoming ray of branch A with the clockwise
    # outgoing ray of branch B, and vice versa.  We record every pairing
    # of an A-ray with a B-ray whose concatenation is anticlockwise.
    terms = []
    a_rays = [r for r in rays if r["branch"] == "A"]
    b_rays = [r for r in rays if r["branch"] == "B"]
    for ra in a_rays:
        for rb in b_rays:
            path = list(reversed(ra["path"])) + rb["path"]
            wind = winding_about(path, puncture)
            # A_f^+ contains positive-length anticlockwise paths only.
            # Winding zero is an idempotent contribution in Mor^x, not a
            # bounding-cochain term; clockwise paths are discarded as well.
            if wind > 0:
                terms.append(("A+B", wind, ra["end"], rb["end"]))
    return {"carriers": 2, "terms": terms, "rays": rays, "d0": 0.0, "d1": 0.0}


def chart_of_orbit(orbs, target):
    hit = locate_orbit(orbs, target)
    if hit is None:
        return None
    _, pp, pre = hit
    xy = pillow_chart(pp[0], pp[1])
    if xy is None:
        return {"xy": None, "face": None, "pre": pre, "P": pp}
    return {"xy": xy, "face": face_of(xy[0], xy[1]), "pre": pre, "P": pp}


def _neighbor_skeleton_hits(hits, param, period, uniq=None, active=None):
    eligible = hits
    if active is not None and uniq is not None:
        eligible = [h for h in hits if _match_unique(h, uniq) in active]
    prev = min(eligible, key=lambda h: (param - h["param"]) % period
               if (param - h["param"]) % period > 1e-9 else period)
    nxt = min(eligible, key=lambda h: (h["param"] - param) % period
              if (h["param"] - param) % period > 1e-9 else period)
    return prev, nxt


def exact_rays_for_preimage(blue, crossing, hits, uniq, active=None,
                            subdivisions=6):
    """The four node-to-skeleton rays at one T^2 lift of a physical node."""
    period = len(blue) - 1
    rays = []
    for branch, k, t in (("A", crossing["kA"], crossing["tA"]),
                         ("B", crossing["kB"], crossing["tB"])):
        param = k + t
        prev, nxt = _neighbor_skeleton_hits(
            hits, param, period, uniq=uniq, active=active)
        backward = list(reversed(
            torus_chart_path_between(
                blue, prev["param"], param, subdivisions=subdivisions)))
        forward = torus_chart_path_between(
            blue, param, nxt["param"], subdivisions=subdivisions)
        for direction, endpoint, path in (("backward", prev, backward),
                                           ("forward", nxt, forward)):
            dot = _match_unique(endpoint, uniq)
            rays.append({
                "branch": branch,
                "direction": direction,
                "param": param,
                "endpoint": endpoint,
                "dot": dot,
                "path": path,
            })
    return rays


def _right_turn(incoming, outgoing, tol=1e-7):
    """Whether an incoming tangent turns clockwise onto an outgoing one."""
    if len(incoming) < 2 or len(outgoing) < 2:
        return False, 0.0
    vin = (incoming[-1][0] - incoming[-2][0],
           incoming[-1][1] - incoming[-2][1])
    vout = (outgoing[1][0] - outgoing[0][0],
            outgoing[1][1] - outgoing[0][1])
    ang = oriented_angle(vin, vout)
    return ang < -tol, ang


def resolve_preimage_exact(blue, crossing, source_branch, target_branch,
                           face, hits, uniq, active=None, subdivisions=6):
    """Resolve one lifted ordered branch jump by exact right-turn paths."""
    rays = exact_rays_for_preimage(
        blue, crossing, hits, uniq, active=active,
        subdivisions=subdivisions)
    source = [r for r in rays if r["branch"] == source_branch]
    target = [r for r in rays if r["branch"] == target_branch]
    candidates = []
    for rs in source:
        incoming = list(reversed(rs["path"]))
        for rt in target:
            right, angle = _right_turn(incoming, rt["path"])
            if not right:
                continue
            path = incoming + rt["path"][1:]
            word = face_path_word(
                path, rs["endpoint"]["arc"], rt["endpoint"]["arc"], face)
            candidates.append({
                "src": rs["dot"],
                "tgt": rt["dot"],
                "word": word,
                "angle": angle,
                "src_ray": (rs["branch"], rs["direction"]),
                "tgt_ray": (rt["branch"], rt["direction"]),
                "src_pt": rs["endpoint"]["pt"],
                "tgt_pt": rt["endpoint"]["pt"],
            })
    positive = [(c["src"], c["word"], c["tgt"])
                for c in candidates if c["word"][1] > 0]
    return {
        "rays": rays,
        "candidates": candidates,
        "u": collect_f2(positive),
    }


def ordered_degree_one(crossing):
    deg = akaho_degrees(crossing)
    return ("A", "B") if deg["A->B"] == 1 else ("B", "A")


def _join_record(pre, endpoints):
    """Positive type-D join with the given unordered endpoint pair."""
    pair = frozenset(endpoints)
    matches = [(i, word, j) for i, word, j in pre["delta"]
               if frozenset((i, j)) == pair]
    if len(matches) != 1:
        return None
    i, word, j = matches[0]
    return {"src": i, "tgt": j, "face": word[0], "length": word[1]}


def physical_carrier_pair(blue, crossing, data):
    """Carrier joins for the ordered degree-one jump at one T^2 lift."""
    pre = data["precurve"]
    active = set(pre["vertices"])
    rays = exact_rays_for_preimage(
        blue, crossing, data["hits"], data["unique"], active=active)
    carriers = {}
    for branch in ("A", "B"):
        dots = [r["dot"] for r in rays if r["branch"] == branch]
        carriers[branch] = _join_record(pre, dots)
    source, target = ordered_degree_one(crossing)
    return carriers.get(source), carriers.get(target)


def _cross2(u, v):
    return u[0] * v[1] - u[1] * v[0]


def _sub2(u, v):
    return u[0] - v[0], u[1] - v[1]


def _open_segment_intersection(a, b, c, d, tol=1e-8):
    """Transverse interior intersection of two segments in R^2."""
    r, s = _sub2(b, a), _sub2(d, c)
    den = _cross2(r, s)
    if abs(den) < tol:
        return None
    ca = _sub2(c, a)
    t = _cross2(ca, s) / den
    u = _cross2(ca, r) / den
    if not (tol < t < 1.0 - tol and tol < u < 1.0 - tol):
        return None
    return t, u, (a[0] + t * r[0], a[1] + t * r[1])


def standard_pairing_intersections(first, second, eps=0.08):
    """KWZ Def. 5.15 intersections and Def. 5.17 resolutions.

    ``first`` and ``second`` are positive two-sided joins.  The first is put
    in the +epsilon position and the second in the -epsilon position in the
    universal polar coordinates of their common face.  Each returned record
    is one lower intersection together with its full right-turn resolution.
    """
    if first is None or second is None or first["face"] != second["face"]:
        return []
    face = first["face"]
    nf = 1 if face in ("L", "R") else 2
    side = lambda dot: (0 if face in ("L", "R")
                        else (0 if dot["arc"] == "L" else 1))
    # Integer starts s and r.  Their reductions mod n_f are the skeleton sides.
    s0 = side(first["src_info"])
    r0 = side(second["src_info"])
    a, b = first["length"], second["length"]
    first_segments = [
        ((s0 + eps, 0.0), (s0 + eps + a, float(a)), "spiral"),
        ((s0 + eps + a, float(a)), (s0 + eps + a, 0.0), "radial"),
    ]
    second_segments = [
        ((r0 - eps, 0.0), (r0 - eps + b, float(b)), "spiral"),
        ((r0 - eps + b, float(b)), (r0 - eps + b, 0.0), "radial"),
    ]
    intersections = []
    for ia, (p0, p1, ptag) in enumerate(first_segments):
        for ib, (q0, q1, qtag) in enumerate(second_segments):
            for shift in range(-a - b - 2, a + b + 3):
                qq0 = (q0[0] + shift * nf, q0[1])
                qq1 = (q1[0] + shift * nf, q1[1])
                hit = _open_segment_intersection(p0, p1, qq0, qq1)
                if hit is None:
                    continue
                _, _, point = hit
                vf, vs = _sub2(p1, p0), _sub2(qq1, qq0)
                first_ends = [
                    (first["src"], s0 + eps, vf, "start"),
                    (first["tgt"], s0 + eps + a, (-vf[0], -vf[1]), "end"),
                ]
                second_ends = [
                    (second["tgt"], r0 - eps + b + shift * nf, vs, "end"),
                    (second["src"], r0 - eps + shift * nf,
                     (-vs[0], -vs[1]), "start"),
                ]
                terms = []
                ignored = []
                for src, sx, vin, stag in first_ends:
                    for tgt, tx, vout, ttag in second_ends:
                        angle = math.atan2(
                            _cross2(vin, vout),
                            vin[0] * vout[0] + vin[1] * vout[1])
                        if angle >= -1e-8:
                            continue
                        length = int(round(tx - sx))
                        rec = (src, (face, length), tgt)
                        if length > 0:
                            terms.append(rec)
                        else:
                            ignored.append(rec)
                intersections.append({
                    "segments": (ptag, qtag),
                    "shift": shift,
                    "point": point,
                    "resolution": collect_f2(terms),
                    "clockwise_or_zero": ignored,
                })
    intersections.sort(key=lambda z: (z["point"][1], z["point"][0]))
    return intersections


def standard_pairing_for_node(blue, preimages, data):
    """Standard carrier-join calculation for one physical node orbit."""
    pairs = [physical_carrier_pair(blue, crossing, data)
             for crossing in preimages]
    sig = lambda j: (None if j is None else
                     (j["src"], j["face"], j["length"], j["tgt"]))
    signatures = [(sig(a), sig(b)) for a, b in pairs]
    if not signatures or any(z != signatures[0] for z in signatures[1:]):
        return {"agree": False, "pairs": pairs, "intersections": []}
    first, second = pairs[0]
    if first is None or second is None:
        return {"agree": False, "pairs": pairs, "intersections": []}
    for join in (first, second):
        join["src_info"] = data["unique"][join["src"]]
        join["tgt_info"] = data["unique"][join["tgt"]]
    return {
        "agree": True,
        "pairs": pairs,
        "first": first,
        "second": second,
        "intersections": standard_pairing_intersections(first, second),
    }


def resolve_orbit_exact(blue, preimages, face, data, subdivisions=6):
    """Apply the right-turn proxy to both deck lifts and require agreement.

    This becomes the KWZ Definition 5.17 resolution only after an explicit
    comparison with first/second pairing position; that comparison is open.
    """
    lifts = []
    active = None
    if "precurve" in data:
        active = set(data["precurve"]["vertices"])
    for crossing in preimages:
        source, target = ordered_degree_one(crossing)
        res = resolve_preimage_exact(
            blue, crossing, source, target, face, data["hits"], data["unique"],
            active=active, subdivisions=subdivisions)
        res["order"] = source + "->" + target
        lifts.append(res)
    lift_terms = [sorted(res["u"], key=repr) for res in lifts]
    agree = bool(lift_terms) and all(ts == lift_terms[0] for ts in lift_terms[1:])
    return {
        "lifts": lifts,
        "agree": agree,
        "u": lift_terms[0] if agree else [],
    }


def report_mc(blue, orbs, data):
    """Resolve all degree-one physical supports on the cancelled precurve."""
    print("== right-turn proxy on the fully cancelled precurve ==")
    pre = data.get("precurve")
    if pre is None:
        check(False, "single-copy precurve is available")
        return None
    delta = pre["delta"]
    active = set(pre["vertices"])
    print(f"  |V|={len(active)}  |δ|={len(delta)}")
    results = {}
    for name, tgt in TARGETS.items():
        ch = chart_of_orbit(orbs, tgt)
        print(f"  -- {name} --")
        if ch is None or ch["xy"] is None or ch["face"] is None:
            print(f"    not in a finite chart face")
            check(False, f"{name} sits in a finite face")
            continue
        face = ch["face"]
        print(f"    chart {tuple(round(v, 4) for v in ch['xy'])}  face={face}")
        mor = resolve_orbit_exact(blue, ch["pre"], face, data)
        def resolution_signature(res):
            return (
                res["agree"],
                tuple(res["u"]),
                tuple(tuple((c["src"], c["word"], c["tgt"])
                            for c in lift["candidates"])
                      for lift in res["lifts"]),
            )
        base_sig = resolution_signature(mor)
        stable = all(
            resolution_signature(resolve_orbit_exact(
                blue, ch["pre"], face, data, subdivisions=ns)) == base_sig
            for ns in (3, 12, 24))
        for j, lift in enumerate(mor["lifts"]):
            cands = [(c["src"], c["word"], c["tgt"])
                     for c in lift["candidates"]]
            print(f"    lift {j} order {lift['order']}: right turns {cands}")
        check(mor["agree"], f"{name} deck lifts give the same proxy morphism")
        check(stable, f"{name} resolution is stable at path subdivisions 3,6,12,24")
        u = mor["u"]
        print(f"    φ in Mor^+: {u if u else 'zero'}")
        check(all(i in active and j in active for i, _, j in u),
              f"{name} morphism uses surviving generators")
        residue, du, uu = mc_residue(delta, u)
        print(f"    |u|={len(u)}  |u²|={len(uu)}  |D(u)|={len(du)}  "
              f"|D(u)+u²|={len(residue)}")
        if residue:
            print(f"    residue: {residue[:8]}"
                  + (" ..." if len(residue) > 8 else ""))
        ok = len(u) > 0 and len(residue) == 0
        print(f"    verdict: {'MC solution' if ok else 'not an MC solution'}")
        results[name] = dict(face=face, u=u, residue=residue, ok=ok)

    grades, grading_ok = generator_bigradings(pre["vertices"], delta)
    check(grading_ok and len(grades) == len(active),
          "the 30-dot differential has a consistent relative KWZ bigrading")
    proxy_degrees = {
        name: [morphism_term_bigrading(term, grades) for term in rec["u"]]
        for name, rec in results.items()
    }
    print(f"  proxy KWZ (q,delta)-degrees: {proxy_degrees}")
    check(proxy_degrees["S25"] == [(Fraction(-8), Fraction(0))]
          and proxy_degrees["S74"] == [(Fraction(4), Fraction(0))],
          "the nonzero proxies have distinct KWZ quantum degrees and delta-degree 0")
    print("  grading warning: the proxy terms are not of differential bidegree "
          "(q,delta)=(0,-1); pairing-position/gradation comparison remains open")

    # Search all distinct nonzero sums of the four proxy morphisms.  Zero
    # proxies do not enlarge the search space.
    names = list(results)
    seen, searched, solutions = set(), [], []
    for mask in range(1, 1 << len(names)):
        support = [names[k] for k in range(len(names)) if mask & (1 << k)]
        u = collect_f2([
            term for name in support for term in results[name]["u"]])
        if not u:
            continue
        key = tuple(sorted(u, key=repr))
        if key in seen:
            continue
        seen.add(key)
        residue, du, uu = mc_residue(delta, u)
        searched.append((support, u, residue))
        if not residue:
            solutions.append((support, u))
    print("  -- small-sum search --")
    for support, u, residue in searched:
        print(f"    {support}: u={u}, residue={residue}")
    check(len(searched) == 3,
          "the two nonzero proxies generate three distinct nonzero sums")
    check(not solutions,
          "no nonzero sum of the four proxy morphisms is Maurer-Cartan")
    results["small_sums"] = dict(searched=searched, solutions=solutions)
    return results


def report_standard_pairing(blue, orbs, data):
    """Run Defs. 5.15 and 5.17 on the four physical carrier pairs."""
    print("== KWZ first/second pairing position on the physical carrier joins ==")
    pre = data["precurve"]
    grades, grading_ok = generator_bigradings(pre["vertices"], pre["delta"])
    check(grading_ok, "relative KWZ bigrading is available for pairing resolutions")
    counts = {}
    all_mc = True
    # Explicit degree-zero primitives, found in the homogeneous End complex.
    # The intersection ordering is radial/spiral first, spiral/radial second.
    primitives = {
        ("S69", 0): [(27, ("1", 0), 20)],
        ("S69", 1): [(27, ("1", 0), 0)],
        ("S18", 0): [(22, ("1", 0), 2), (23, ("1", 0), 3)],
        ("S18", 1): [(24, ("1", 0), 3), (27, ("M", 1), 2)],
    }
    all_gauge_trivial = True
    standard_cycles = []
    for name, tgt in TARGETS.items():
        hit = locate_orbit(orbs, tgt)
        print(f"  -- {name} --")
        if hit is None:
            check(False, f"{name} is available for the standard-pairing calculation")
            continue
        _, _, crossings = hit
        result = standard_pairing_for_node(blue, crossings, data)
        check(result["agree"], f"{name} deck lifts determine the same ordered carrier pair")
        if not result["agree"]:
            continue
        first, second = result["first"], result["second"]
        fmt = lambda j: f"{j['src']} --{j['face']}^{j['length']}--> {j['tgt']}"
        print(f"    degree-one carrier map: {fmt(first)}  to  {fmt(second)}")
        xs = result["intersections"]
        pairing_sigs = []
        for epsilon in (0.02, 0.05, 0.08, 0.12):
            trial = standard_pairing_intersections(first, second, eps=epsilon)
            pairing_sigs.append(tuple(
                (x["segments"], tuple(sorted(x["resolution"], key=repr)))
                for x in trial))
        check(all(sig == pairing_sigs[0] for sig in pairing_sigs[1:]),
              f"{name} standard resolution is stable for four epsilon choices")
        counts[name] = len(xs)
        print(f"    lower intersections in standard position: {len(xs)}")
        for k, x in enumerate(xs):
            u = x["resolution"]
            degrees = [morphism_term_bigrading(term, grades) for term in u]
            residue, du, uu = mc_residue(pre["delta"], u)
            all_mc = all_mc and bool(u) and not residue
            print(f"      x{k}: {x['segments']}  phi={u}  "
                  f"degrees={degrees}  |Dphi|={len(du)}  "
                  f"|phi^2|={len(uu)}  |MC|={len(residue)}")
            check(bool(u), f"{name} standard intersection x{k} has nonzero resolution")
            check(len(set(degrees)) <= 1,
                  f"{name} standard intersection x{k} resolves homogeneously")
            check(not residue,
                  f"{name} standard intersection x{k} is Maurer-Cartan algebraically")
            v = primitives.get((name, k))
            if v is not None:
                dv = collect_f2(
                    apply_mul(pre["delta"], v) + apply_mul(v, pre["delta"]))
                vv, vu, uv = apply_mul(v, v), apply_mul(v, u), apply_mul(u, v)
                gauge_residue = collect_f2(dv + u + vu)
                trivial = (dv == u and not vv and not vu and not uv
                           and not gauge_residue)
                all_gauge_trivial = all_gauge_trivial and trivial
                print(f"        primitive v={v}  Dv={dv}  "
                      f"|v^2|={len(vv)} |vu|={len(vu)} |uv|={len(uv)}")
                check(trivial,
                      f"{name} standard intersection x{k} is gauge-trivial via 1+v")
                standard_cycles.append((name, k, u, v))
    check(counts == {"S69": 2, "S18": 2, "S25": 0, "S74": 0},
          "standard pairing gives two R-face intersections and no M(2)-to-M(1) intersections")
    check(all_mc, "all four nonzero standard resolutions are D-closed and square-zero")
    check(all_gauge_trivial,
          "all four Maurer-Cartan deformations are isomorphic to the undeformed type-D object")
    all_sum_trivial = True
    for mask in range(1, 1 << len(standard_cycles)):
        u = collect_f2([
            term for k, (_, _, cycle, _) in enumerate(standard_cycles)
            if mask & (1 << k) for term in cycle])
        v = collect_f2([
            term for k, (_, _, _, primitive) in enumerate(standard_cycles)
            if mask & (1 << k) for term in primitive])
        dv = collect_f2(
            apply_mul(pre["delta"], v) + apply_mul(v, pre["delta"]))
        all_sum_trivial = all_sum_trivial and (
            not collect_f2(dv + u)
            and not apply_mul(v, v)
            and not apply_mul(v, u)
            and not apply_mul(u, v))
    check(all_sum_trivial,
          "every nonzero sum of the four standard cycles is gauge-trivial")
    print("  continuation warning: pairing position does not by itself identify either "
          "standard R-face intersection with the original physical node")
    print("  grading warning: the standard resolutions have delta-degree -1 but "
          "quantum degree 6 or 8, not the quantum degree of the type-D differential")
    return counts


def report_end_cohomology(data):
    """Degree-minus-one cohomology of the corrected blue End complex."""
    print("== graded degree-minus-one cohomology of End(N) ==")
    pre = data["precurve"]
    grades, grading_ok = generator_bigradings(pre["vertices"], pre["delta"])
    check(grading_ok, "relative grading is consistent for the End calculation")
    bases = {
        degree: end_terms_in_delta_degree(data, pre, grades, degree)
        for degree in (0, -1, -2)
    }
    by_q = {degree: {} for degree in bases}
    for degree, basis in bases.items():
        for term in basis:
            qdegree = morphism_term_bigrading(term, grades)[0]
            by_q[degree].setdefault(qdegree, []).append(term)
    rows = {}
    qvalues = sorted(by_q[-1])
    all_closed = True
    for qdegree in qvalues:
        b0 = by_q[0].get(qdegree, [])
        b1 = by_q[-1].get(qdegree, [])
        b2 = by_q[-2].get(qdegree, [])
        ix1 = {term: k for k, term in enumerate(b1)}
        ix2 = {term: k for k, term in enumerate(b2)}
        cols0, cols1 = [], []
        for term in b0:
            image = end_differential(pre["delta"], [term])
            all_closed = all_closed and all(z in ix1 for z in image)
            cols0.append([ix1[z] for z in image])
        for term in b1:
            image = end_differential(pre["delta"], [term])
            all_closed = all_closed and all(z in ix2 for z in image)
            cols1.append([ix2[z] for z in image])
        rank0, rank1 = sparse_column_rank(cols0), sparse_column_rank(cols1)
        rows[qdegree] = len(b1) - rank0 - rank1
    check(all_closed, "the fixed-delta bases are closed under the End differential")
    nonzero = {int(q): dim for q, dim in rows.items() if dim}
    print(f"  nonzero H^(q,-1) dimensions: {nonzero}")
    expected = {
        -16: 1, -14: 1, -12: 1, -10: 1, -8: 2, -6: 1,
        -2: 5, 0: 1, 2: 2, 4: 1, 6: 2, 8: 1, 10: 2, 12: 1,
    }
    check(nonzero == expected,
          "the degree-minus-one End cohomology profile matches the exact sparse calculation")
    representatives = end_h_representatives()
    reps_ok = {int(q): len(reps) for q, reps in representatives.items()} == expected
    for qdegree, reps in representatives.items():
        target = by_q[-1][qdegree]
        target_index = {term: k for k, term in enumerate(target)}
        boundary_columns = []
        for term in by_q[0].get(qdegree, []):
            image = end_differential(pre["delta"], [term])
            boundary_columns.append([target_index[z] for z in image])
        rank_boundary = sparse_column_rank(boundary_columns)
        rep_columns = []
        for rep in reps:
            reps_ok = reps_ok and not end_differential(pre["delta"], rep)
            rep_columns.append([target_index[z] for z in rep])
        reps_ok = reps_ok and sparse_column_rank(
            boundary_columns + rep_columns) == rank_boundary + len(reps)
    check(reps_ok,
          "the 22 displayed homogeneous cycles form a cohomology basis")

    # Quantum degree zero has one class.  Its short representative is the
    # horizontal resolution between the M^2 joins 7->4 and 12->15.
    qzero0 = by_q[0].get(Fraction(0), [])
    qzero1 = by_q[-1].get(Fraction(0), [])
    ix = {term: k for k, term in enumerate(qzero1)}
    columns = []
    for term in qzero0:
        image = end_differential(pre["delta"], [term])
        columns.append([ix[z] for z in image])
    u0 = [(7, ("M", 2), 12), (4, ("M", 2), 15)]
    u0_column = [ix[z] for z in u0]
    rank_before = sparse_column_rank(columns)
    rank_after = sparse_column_rank(columns + [u0_column])
    check(not end_differential(pre["delta"], u0) and not apply_mul(u0, u0),
          "u0=(7 --M^2--> 12)+(4 --M^2--> 15) is closed and square-zero")
    check(rows[Fraction(0)] == 1 and rank_after == rank_before + 1,
          "[u0] spans the only quantum-preserving degree-minus-one End class")
    delta_primitive = [
        (i, ("1", 0), i)
        for i in (19, 17, 15, 11, 9, 7, 3, 1,
                  20, 22, 24, 28, 30, 32, 34)
    ]
    check(not collect_f2(end_differential(pre["delta"], delta_primitive)
                         + pre["delta"]),
          "delta itself is exact via the alternating idempotent projector")
    print("  unique q=0 MC class: u0=[(7,M^2,12),(4,M^2,15)] "
          "(physical carrier orbit S23)")

    # The unique cohomology representative found below with the desired
    # red-pairing rank direction.  It is a wrapped standard-position class,
    # not the image of a physical self-crossing of the original PL curve.
    ustar = [(1, ("M", 1), 27), (0, ("M", 2), 24)]
    uaux = [(1, ("M", 1), 33)]
    ucombo = collect_f2(ustar + uaux)
    ustar_degree = [morphism_term_bigrading(z, grades) for z in ustar]
    check(ustar_degree == [(Fraction(-8), Fraction(-1))] * 2
          and not end_differential(pre["delta"], ustar)
          and not apply_mul(ustar, ustar),
          "u* is a homogeneous degree-(-8,-1), closed, square-zero End cycle")
    qm8_source = by_q[0][Fraction(-8)]
    qm8_target = by_q[-1][Fraction(-8)]
    qm8_ix = {term: k for k, term in enumerate(qm8_target)}
    qm8_columns = [
        [qm8_ix[z] for z in end_differential(pre["delta"], [term])]
        for term in qm8_source
    ]
    qm8_rank = sparse_column_rank(qm8_columns)
    qm8_with_ustar = sparse_column_rank(
        qm8_columns + [[qm8_ix[z] for z in ustar]])
    check(qm8_with_ustar == qm8_rank + 1,
          "u* is not an End boundary")
    first = {"src": 1, "tgt": 0, "face": "M", "length": 1,
             "src_info": data["unique"][1], "tgt_info": data["unique"][0]}
    second = {"src": 27, "tgt": 24, "face": "M", "length": 2,
              "src_info": data["unique"][27], "tgt_info": data["unique"][24]}
    wrapped = standard_pairing_intersections(first, second)
    check(any(x["resolution"] == ustar for x in wrapped),
          "u* is the horizontal resolution of the M^1-to-M^2 wrapped intersection")
    check(not end_differential(pre["delta"], uaux)
          and not apply_mul(uaux, uaux)
          and not apply_mul(ucombo, ucombo),
          "the second q=-8 class and its sum with u* are also MC cycles")
    print("  +2 Hom-shift coset representative: "
          "u*=[(1,M,27),(0,M^2,24)]; neither it nor u*+(1,M,33) "
          "has a physical carrier crossing")
    return rows


def _q7_red_type_d(red):
    """The six-dot red cycle after its unique skeleton-bigon cancellation."""
    data = encode_type_d(red)
    check(len(data["hits"]) == len(data["unique"]) == 8,
          "the q=7 red earring has eight skeleton dots before cancellation")
    pre = cyclic_precurve(data)
    check(len(pre["vertices"]) == 6 and len(pre["delta"]) == 6
          and len(pre["cancellations"]) == 1,
          "the red six-cycle is derived by one zero-join cancellation")
    check(not pre["residue"],
          "the six-dot red type-D differential squares to zero")
    return data, pre["vertices"], pre["delta"]


def _bounded_hom_basis(red_data, red_vertices, blue_data, blue_pre, cutoff):
    """All red-to-blue Hom terms whose face-word length is at most cutoff."""
    terms = []
    for i in red_vertices:
        arc_i = red_data["unique"][i]["arc"]
        side_i = 0 if arc_i == "L" else 1
        for j in blue_pre["vertices"]:
            arc_j = blue_data["unique"][j]["arc"]
            side_j = 0 if arc_j == "L" else 1
            if arc_i == arc_j:
                terms.append((i, ("1", 0), j))
            if arc_i == arc_j == "L":
                terms.extend((i, ("L", m), j) for m in range(1, cutoff + 1))
            if arc_i == arc_j == "R":
                terms.extend((i, ("R", m), j) for m in range(1, cutoff + 1))
            terms.extend(
                (i, ("M", m), j) for m in range(1, cutoff + 1)
                if (side_i + m - side_j) % 2 == 0)
    return terms


def _hom_differential_column(term, red_delta, blue_delta):
    """The exact type-D Hom differential of one red-to-blue morphism."""
    i, word, j = term
    image = []
    for src, left_word, tgt in red_delta:
        if tgt != i:
            continue
        product = mul_word(left_word, word)
        if product is not None:
            image.append((src, product, j))
    for src, right_word, tgt in blue_delta:
        if src != j:
            continue
        product = mul_word(word, right_word)
        if product is not None:
            image.append((i, product, tgt))
    return collect_f2(image)


def _bounded_hom_matrix(red_data, red_vertices, red_delta,
                        blue_data, blue_pre, blue_delta, cutoff):
    """Exact D on inputs of length <= cutoff, with every output retained."""
    domain = _bounded_hom_basis(
        red_data, red_vertices, blue_data, blue_pre, cutoff)
    max_step = max(word[1] for _, word, _ in red_delta + blue_delta)
    target = _bounded_hom_basis(
        red_data, red_vertices, blue_data, blue_pre, cutoff + max_step)
    index = {term: k for k, term in enumerate(target)}
    columns = []
    for term in domain:
        image = _hom_differential_column(term, red_delta, blue_delta)
        check_target = [z for z in image if z not in index]
        if check_target:
            raise AssertionError(f"Hom output escaped exact target: {check_target}")
        columns.append([index[z] for z in image])
    return domain, target, columns


def finite_support_homology(red_data, red_vertices, red_delta,
                            blue_data, blue_pre, blue_delta, cutoff,
                            primitive_extra=0):
    """Exact low-support cycles modulo boundaries from a bounded primitive set.

    Cycles are elements supported in word lengths <= ``cutoff`` whose full
    differential vanishes; outputs above the cutoff are therefore not silently
    discarded.  Boundary dimension is computed as

        dim(im D|C_<=K intersect C_<=cutoff)
          = rank(D|C_<=K) - rank(high-output projection),

    where K=cutoff+primitive_extra.  Increasing K tests whether cancellations
    among longer primitives create any additional low-support boundary.
    """
    cycle_domain, _, cycle_columns = _bounded_hom_matrix(
        red_data, red_vertices, red_delta,
        blue_data, blue_pre, blue_delta, cutoff)
    cycle_rank = sparse_column_rank(cycle_columns)
    cycle_dim = len(cycle_domain) - cycle_rank

    primitive_cutoff = cutoff + primitive_extra
    _, boundary_target, boundary_columns = _bounded_hom_matrix(
        red_data, red_vertices, red_delta,
        blue_data, blue_pre, blue_delta, primitive_cutoff)
    full_rank = sparse_column_rank(boundary_columns)
    high_rows = {
        k for k, (_, word, _) in enumerate(boundary_target)
        if word[1] > cutoff
    }
    high_rank = sparse_column_rank(
        [[row for row in column if row in high_rows]
         for column in boundary_columns])
    boundary_dim = full_rank - high_rank
    homology_dim = cycle_dim - boundary_dim
    if homology_dim < 0:
        raise AssertionError("boundary intersection exceeds the cycle space")
    return {
        "chain": len(cycle_domain),
        "cycles": cycle_dim,
        "boundaries": boundary_dim,
        "homology": homology_dim,
    }


def _poly_mul_f2(a, b):
    """Multiply bit-encoded polynomials in F2[U]."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
    return out


def _tower_base_length(face, side_i, side_j):
    if face in ("L", "R"):
        return 1
    return 1 if (side_i + 1 - side_j) % 2 == 0 else 2


def _face_tower_complex(face, red_data, red_vertices, red_delta,
                        blue_data, blue_pre, blue_delta):
    """Positive face-word Hom complex as a matrix over F2[U].

    For L and R, multiplication by U raises word length by one.  For M it
    raises length by two, preserving the endpoint parity condition.
    """
    period = 1 if face in ("L", "R") else 2
    labels = []
    for i in red_vertices:
        arc_i = red_data["unique"][i]["arc"]
        side_i = 0 if arc_i == "L" else 1
        for j in blue_pre["vertices"]:
            arc_j = blue_data["unique"][j]["arc"]
            side_j = 0 if arc_j == "L" else 1
            if face == "L" and not arc_i == arc_j == "L":
                continue
            if face == "R" and not arc_i == arc_j == "R":
                continue
            labels.append(
                (i, j, _tower_base_length(face, side_i, side_j)))
    index = {(i, j): k for k, (i, j, _) in enumerate(labels)}
    matrix = {}
    for col, (i, j, base) in enumerate(labels):
        images = []
        for src, word, tgt in red_delta:
            if tgt == i and word[0] == face:
                images.append((src, j, base + word[1]))
        for src, word, tgt in blue_delta:
            if src == j and word[0] == face:
                images.append((i, tgt, base + word[1]))
        for out_i, out_j, length in images:
            row = index[(out_i, out_j)]
            out_base = labels[row][2]
            if length < out_base or (length - out_base) % period:
                raise AssertionError("face-tower endpoint parity mismatch")
            coeff = 1 << ((length - out_base) // period)
            key = (row, col)
            matrix[key] = matrix.get(key, 0) ^ coeff
            if not matrix[key]:
                del matrix[key]
    return labels, matrix


def _identity_tail_vectors(face, labels, identities, red_delta, blue_delta):
    """Images of identity morphisms in one positive face-word tower."""
    period = 1 if face in ("L", "R") else 2
    index = {(i, j): (k, base) for k, (i, j, base) in enumerate(labels)}
    vectors = []
    for term in identities:
        vector = {}
        for i, word, j in _hom_differential_column(
                term, red_delta, blue_delta):
            if word[0] != face:
                continue
            row, base = index[(i, j)]
            if word[1] < base or (word[1] - base) % period:
                raise AssertionError("identity image misses its face tower")
            coeff = 1 << ((word[1] - base) // period)
            vector[row] = vector.get(row, 0) ^ coeff
            if not vector[row]:
                del vector[row]
        vectors.append(vector)
    return vectors


def _polynomial_image(matrix, vector):
    out = {}
    for (row, col), coeff in matrix.items():
        if col not in vector:
            continue
        out[row] = out.get(row, 0) ^ _poly_mul_f2(
            coeff, vector[col])
        if not out[row]:
            del out[row]
    return out


def _cancel_unit_tower_pairs(labels, matrix, beta_vectors):
    """Cancel every unit arrow and transport the identity-image cycles."""
    active = set(range(len(labels)))
    cancelled = []
    while True:
        pivot = next(((row, col) for (row, col), coeff in matrix.items()
                      if coeff == 1 and row in active and col in active
                      and row != col), None)
        if pivot is None:
            break
        pivot_row, pivot_col = pivot
        col_entries = [
            (row, coeff) for (row, col), coeff in matrix.items()
            if col == pivot_col and row in active and row != pivot_row
        ]
        row_entries = [
            (col, coeff) for (row, col), coeff in matrix.items()
            if row == pivot_row and col in active and col != pivot_col
        ]
        for row, left in col_entries:
            for col, right in row_entries:
                key = (row, col)
                matrix[key] = matrix.get(key, 0) ^ _poly_mul_f2(
                    left, right)
                if not matrix[key]:
                    del matrix[key]
        # Project each cycle away from the contractible pivot pair.  First
        # apply the source-basis changes e_col' = e_col + right*e_pivot_col
        # which remove the other entries in the pivot row.  In the new basis
        # the pivot-column coordinate therefore gains right times the old
        # e_col coordinate.  Then adding the pivot-row coefficient times
        # D(pivot_col) removes the pivot-row component.
        for vector in beta_vectors:
            source_coeff = vector.get(pivot_col, 0)
            for col, right in row_entries:
                if col in vector:
                    source_coeff ^= _poly_mul_f2(right, vector[col])
            if source_coeff:
                vector[pivot_col] = source_coeff
            else:
                vector.pop(pivot_col, None)
            pivot_coeff = vector.get(pivot_row, 0)
            if pivot_coeff:
                for row, coeff in col_entries:
                    vector[row] = vector.get(row, 0) ^ _poly_mul_f2(
                        coeff, pivot_coeff)
                    if not vector[row]:
                        vector.pop(row, None)
            vector.pop(pivot_row, None)
            if vector.get(pivot_col, 0):
                raise AssertionError(
                    "identity-image cycle retains a cancelled source")
            vector.pop(pivot_col, None)
        active.remove(pivot_row)
        active.remove(pivot_col)
        for key in [key for key in matrix
                    if key[0] in pivot or key[1] in pivot]:
            del matrix[key]
        cancelled.append(pivot)
    return active, matrix, beta_vectors, cancelled


def exact_red_blue_homology(red_data, red_vertices, red_delta,
                            blue_data, blue_pre, blue_delta):
    """Full Hom homology via finite F2[U] tail reduction.

    The positive-word complex is a finite free F2[U]-module in each face.
    Unit cancellation reduces it here to ``U*A`` with ``A^2=0`` and
    ``rank(A)=dim(A)/2``.  Its homology is consequently coker(U)^{rank(A)}.
    The full Hom complex is the cone of the identity-component map into that
    tail homology, so its dimension is finite and computed without a cutoff.
    """
    identities = _bounded_hom_basis(
        red_data, red_vertices, blue_data, blue_pre, 0)
    beta_columns = [[] for _ in identities]
    tail_homology = 0
    ambient_offset = 0
    face_certificate = {}
    for face in ("L", "M", "R"):
        labels, matrix = _face_tower_complex(
            face, red_data, red_vertices, red_delta,
            blue_data, blue_pre, blue_delta)
        vectors = _identity_tail_vectors(
            face, labels, identities, red_delta, blue_delta)
        active, residual, vectors, cancelled = _cancel_unit_tower_pairs(
            labels, dict(matrix), vectors)
        if any(_polynomial_image(residual, vector) for vector in vectors):
            raise AssertionError("transported identity image is not a tail cycle")
        # The residual differential must be exactly U times a constant matrix.
        if any(coeff != 2 for coeff in residual.values()):
            raise AssertionError("tail did not reduce to a U-linear differential")
        order = {old: new for new, old in enumerate(sorted(active))}
        columns = [[] for _ in active]
        for (row, col), coeff in residual.items():
            if coeff == 2:
                columns[order[col]].append(order[row])
        constant_rank = sparse_column_rank(columns)
        if 2 * constant_rank != len(active):
            raise AssertionError("tail has a nonzero free F2[U] homology summand")
        # residual^2=0 is equivalent to A^2=0 after factoring out U.
        a_squared = []
        for column in columns:
            image = 0
            for middle in column:
                for row in columns[middle]:
                    image ^= 1 << row
            a_squared.append(image)
        if any(a_squared):
            raise AssertionError("residual constant matrix does not square to zero")
        for col, vector in enumerate(vectors):
            # In H(U*A), every positive U coefficient is a boundary because
            # ker(A)=im(A); only the constant coefficient remains.
            beta_columns[col].extend(
                ambient_offset + order[row]
                for row, polynomial in vector.items() if polynomial & 1)
        face_certificate[face] = {
            "tower_rank": len(labels),
            "unit_pairs": len(cancelled),
            "residual_rank": len(active),
            "homology": constant_rank,
        }
        tail_homology += constant_rank
        ambient_offset += len(active)
    beta_rank = sparse_column_rank(beta_columns)
    homology = len(identities) + tail_homology - 2 * beta_rank
    return {
        "identities": len(identities),
        "tail_homology": tail_homology,
        "beta_rank": beta_rank,
        "homology": homology,
        "faces": face_certificate,
    }


def _bit_span_basis(vectors):
    """Row-echelon basis for bit-encoded vectors over F2."""
    pivots = {}
    for vector in vectors:
        value = vector
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return pivots


def _bit_in_span(vector, pivots):
    value = vector
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return False
        value ^= pivots[pivot]
    return True


def _kernel_bits(columns, nrows):
    """Basis of the kernel of a bit-column matrix."""
    equations = []
    for row in range(nrows):
        equation = 0
        for col, value in enumerate(columns):
            if value & (1 << row):
                equation ^= 1 << col
        if equation:
            equations.append(equation)
    pivots = _bit_span_basis(equations)
    free = [i for i in range(len(columns)) if i not in pivots]
    kernel = []
    for free_index in free:
        value = 1 << free_index
        for pivot in sorted(pivots):
            if (pivots[pivot] & value).bit_count() % 2:
                value ^= 1 << pivot
        kernel.append(value)
    return kernel


def end_cohomology_data(data, quantum_degree=Fraction(0)):
    """One fixed-quantum-degree End cohomology calculation and representatives."""
    pre = data["precurve"]
    grades, grading_ok = generator_bigradings(pre["vertices"], pre["delta"])
    if not grading_ok:
        raise AssertionError("inconsistent relative grading")
    bases = {
        degree: [
            term for term in end_terms_in_delta_degree(data, pre, grades, degree)
            if morphism_term_bigrading(term, grades)[0] == quantum_degree
        ]
        for degree in (0, -1, -2)
    }
    index1 = {term: i for i, term in enumerate(bases[-1])}
    index2 = {term: i for i, term in enumerate(bases[-2])}
    columns0 = [
        sum(1 << index1[z] for z in end_differential(pre["delta"], [term]))
        for term in bases[0]
    ]
    columns1 = [
        sum(1 << index2[z] for z in end_differential(pre["delta"], [term]))
        for term in bases[-1]
    ]
    boundaries = _bit_span_basis(columns0)
    kernel = _kernel_bits(columns1, len(bases[-2]))
    representatives = []
    span = dict(boundaries)
    for vector in kernel:
        if _bit_in_span(vector, span):
            continue
        representatives.append([
            bases[-1][i] for i in range(len(bases[-1]))
            if vector & (1 << i)
        ])
        span = _bit_span_basis(list(span.values()) + [vector])
    return {
        "grades": grades,
        "bases": bases,
        "rank0": len(boundaries),
        "rank1": sparse_column_rank(
            [[i for i in range(len(bases[-2])) if value & (1 << i)]
             for value in columns1]),
        "representatives": representatives,
    }


def _match_component_dots(component_data, original_data, tolerance=1e-6):
    """Canonical skeleton-dot matching for a surgery supported off the skeleton."""
    matching = {}
    for i, dot in enumerate(component_data["unique"]):
        candidates = [
            (abs(dot["y"] - other["y"]), j)
            for j, other in enumerate(original_data["unique"])
            if dot["arc"] == other["arc"]
        ]
        if not candidates:
            raise AssertionError("component dot has no skeleton-side match")
        distance, j = min(candidates)
        if distance > tolerance:
            raise AssertionError(
                f"surgery moved skeleton dot {i} by {distance}")
        matching[i] = j
    if len(set(matching.values())) != len(matching):
        raise AssertionError("skeleton-dot matching is not injective")
    return matching


def s69_smoothing_deformation(red, blue, target, data, eps=0.006):
    """Directly encode S69 smoothing on the original skeleton-dot module."""
    from bigons import edges_of
    from surgery_check import (bigon_matrix_multi, iota_pairing,
                               locate_supports, smooth_curve)

    candidates = locate_supports(blue, [target])[0]
    support = min(candidates, key=lambda item: _tdist(item[0], target))
    c1, c2 = support[1]
    p2, mismatch = iota_pairing(edges_of(blue), c1, 1, c2, eps)
    if p2 != 1 or mismatch > 1e-8:
        raise AssertionError("S69 smoothing is not iota-compatible")
    components = smooth_curve(blue, [(c1, 1), (c2, p2)], eps)
    classes = [deck_of_path(component) for component in components]
    if classes != [(17, 8), (-2, -1), (-2, -1)]:
        raise AssertionError(f"unexpected S69 component classes {classes}")

    red_data, red_vertices, red_delta = _q7_red_type_d(red)
    encoded = []
    for component in components:
        component_data = encode_type_d(component)
        component_pre = (
            cyclic_precurve(component_data)
            if len(component_data["hits"]) == len(component_data["unique"])
            else geometric_precurve(component_data)
        )
        encoded.append((component_data, component_pre))
    long_index = classes.index((17, 8))
    short_index = next(
        i for i, (_, pre) in enumerate(encoded)
        if pre["vertices"] == red_vertices and pre["delta"] == red_delta
    )

    mapped_vertices = []
    mapped_delta = []
    for index in (long_index, short_index):
        component_data, component_pre = encoded[index]
        matching = _match_component_dots(component_data, data)
        mapped_vertices.extend(matching[i] for i in component_pre["vertices"])
        mapped_delta.extend(
            (matching[i], word, matching[j])
            for i, word, j in component_pre["delta"]
        )
    mapped_delta = collect_f2(mapped_delta)
    if (len(mapped_vertices) != len(set(mapped_vertices))
            or set(mapped_vertices) != set(data["precurve"]["vertices"])):
        raise AssertionError("smoothed components do not match the original module")
    deformation = collect_f2(data["precurve"]["delta"] + mapped_delta)

    long_data, long_pre = encoded[long_index]
    long_hom = exact_red_blue_homology(
        red_data, red_vertices, red_delta,
        long_data, long_pre, long_pre["delta"])["homology"]
    short_hom = exact_red_blue_homology(
        red_data, red_vertices, red_delta,
        red_data, {"vertices": red_vertices}, red_delta)["homology"]
    finite_generators, finite_delta = bigon_matrix_multi(red, components)
    finite_hom = len(finite_generators) - 2 * rank_f2(finite_delta)
    return {
        "components": components,
        "classes": classes,
        "long": long_hom,
        "short": short_hom,
        "finite": finite_hom,
        "delta": mapped_delta,
        "b": deformation,
    }


def physical_smoothing_pairing(red, blue, support, data, pairing=1, eps=0.006,
                               red_type_d=None):
    """Exact wrapped pairing after one iota-equivariant physical smoothing.

    ``support`` is one pillowcase self-intersection orbit ``(point, lifts)``
    returned by :func:`orbit_group`.  Components on the torus which are
    exchanged by the pillowcase involution project to one component and are
    counted only once.  When the canceled skeleton modules before and after
    surgery agree, the function also returns their Maurer--Cartan difference.
    """
    from collections import defaultdict

    from bigons import edges_of
    from surgery_check import bigon_matrix_multi, iota_pairing, smooth_curve

    point, lifts = support
    if len(lifts) != 2:
        raise AssertionError(
            f"physical orbit at {point} has {len(lifts)} torus lifts")
    p2, mismatch = iota_pairing(edges_of(blue), lifts[0], pairing,
                                lifts[1], eps)
    if mismatch > 1e-8:
        raise AssertionError(
            f"smoothing at {point} is not iota-compatible: {mismatch}")
    components = smooth_curve(
        blue, [(lifts[0], pairing), (lifts[1], p2)], eps)

    encoded = []
    quotient_groups = defaultdict(list)
    for index, component in enumerate(components):
        component_data = encode_type_d(component)
        invariant = len(component_data["hits"]) != len(component_data["unique"])
        component_pre = (
            geometric_precurve(component_data)
            if invariant else cyclic_precurve(component_data)
        )
        signature = tuple(sorted(
            (dot["arc"], round(dot["y"], 7))
            for dot in component_data["unique"]
        ))
        encoded.append({
            "index": index,
            "class": deck_of_path(component),
            "data": component_data,
            "precurve": component_pre,
            "invariant": invariant,
            "signature": signature,
        })
        quotient_groups[(invariant, signature)].append(index)

    selected = []
    for (invariant, _), indices in quotient_groups.items():
        if invariant:
            selected.extend(indices)
        else:
            if len(indices) != 2:
                raise AssertionError(
                    "non-invariant torus components do not form an iota pair")
            selected.append(indices[0])

    red_data, red_vertices, red_delta = (
        _q7_red_type_d(red) if red_type_d is None else red_type_d)
    homology = []
    homology_certificates = []
    mapped_vertices = []
    mapped_delta = []
    for index in selected:
        record = encoded[index]
        component_data = record["data"]
        component_pre = record["precurve"]
        certificate = exact_red_blue_homology(
            red_data, red_vertices, red_delta,
            component_data, component_pre, component_pre["delta"]
        )
        homology.append(certificate["homology"])
        homology_certificates.append(certificate)
        matching = _match_component_dots(component_data, data)
        mapped_vertices.extend(
            matching[i] for i in component_pre["vertices"])
        mapped_delta.extend(
            (matching[i], word, matching[j])
            for i, word, j in component_pre["delta"]
        )

    mapped_delta = collect_f2(mapped_delta)
    module_matches = (
        len(mapped_vertices) == len(set(mapped_vertices))
        and set(mapped_vertices) == set(data["precurve"]["vertices"])
    )
    deformation = None
    if module_matches:
        deformation = collect_f2(data["precurve"]["delta"] + mapped_delta)
        if mc_residue(data["precurve"]["delta"], deformation)[0]:
            raise AssertionError("direct physical smoothing is not Maurer--Cartan")

    finite_generators, finite_delta = bigon_matrix_multi(red, components)
    finite_homology = len(finite_generators) - 2 * rank_f2(finite_delta)
    return {
        "point": point,
        "pairings": (pairing, p2),
        "components": components,
        "classes": [record["class"] for record in encoded],
        "selected": selected,
        "wrapped_summands": homology,
        "wrapped_certificates": homology_certificates,
        "wrapped": sum(homology),
        "finite": finite_homology,
        "module_matches": module_matches,
        "delta": mapped_delta if module_matches else None,
        "b": deformation,
    }


def report_admissible_kwz(red, blue, targets, data):
    """Correct q=7 KWZ calculation using the CHKK special puncture (0,pi)."""
    print("== admissible-special-puncture KWZ pairing and S69 deformation ==")
    red_data, red_vertices, red_delta = _q7_red_type_d(red)
    pre = data["precurve"]
    base = exact_red_blue_homology(
        red_data, red_vertices, red_delta, data, pre, pre["delta"])
    print(f"  undeformed exact wrapped Hom dimension={base['homology']}; "
          f"identity/tail/beta={base['identities']}/"
          f"{base['tail_homology']}/{base['beta_rank']}")
    check(base["homology"] == 7
          and (base["identities"], base["tail_homology"], base["beta_rank"])
          == (104, 105, 101),
          "the admissible-chart undeformed wrapped pairing has dimension 7")

    qzero = end_cohomology_data(data, Fraction(0))
    check(len(qzero["representatives"]) == 1,
          "quantum-preserving degree-minus-one End cohomology is one-dimensional")
    u0 = qzero["representatives"][0]
    u0_result = exact_red_blue_homology(
        red_data, red_vertices, red_delta, data, pre,
        collect_f2(pre["delta"] + u0))
    print(f"  q=0 representative={u0}; deformed Hom={u0_result['homology']}")
    check(u0 == [(5, ("M", 1), 10)]
          and not mc_residue(pre["delta"], u0)[0]
          and u0_result["homology"] == 5,
          "the unique quantum-preserving MC class changes 7 to 5")

    grades = qzero["grades"]
    orbits = orbit_group(blue)
    physical_dimensions = {}
    for name, target in targets.items():
        hit = locate_orbit(orbits, target)
        if hit is None:
            raise AssertionError(f"missing translated physical orbit {name}")
        result = standard_pairing_for_node(blue, hit[2], data)
        dimensions = []
        for intersection in result["intersections"]:
            cycle = intersection["resolution"]
            degrees = {morphism_term_bigrading(term, grades) for term in cycle}
            if len(degrees) != 1 or mc_residue(pre["delta"], cycle)[0]:
                raise AssertionError(f"invalid standard resolution at {name}")
            dimensions.append(exact_red_blue_homology(
                red_data, red_vertices, red_delta, data, pre,
                collect_f2(pre["delta"] + cycle))["homology"])
        physical_dimensions[name] = dimensions
    print(f"  standard physical-resolution Hom dimensions={physical_dimensions}")
    check(physical_dimensions == {
        "S69": [7, 7], "S18": [7, 7],
        "S25": [7, 7], "S74": [7, 7],
    }, "all eight single standard resolutions preserve wrapped rank 7")

    smooth = s69_smoothing_deformation(red, blue, targets["S69"], data)
    expected_b = {
        (0, ("R", 1), 19), (26, ("R", 1), 27),
        (0, ("R", 1), 27), (26, ("R", 1), 19),
    }
    b = smooth["b"]
    b_bigradings = [morphism_term_bigrading(term, grades) for term in b]
    degrees = sorted(degree[0] for degree in b_bigradings)
    deformed = exact_red_blue_homology(
        red_data, red_vertices, red_delta, data, pre, smooth["delta"])
    print(f"  S69 switch b={b}; quantum degrees={degrees}")
    print(f"  S69 components: classes={smooth['classes']}; "
          f"wrapped={smooth['long']}+{smooth['short']}={deformed['homology']}; "
          f"finite={smooth['finite']}")
    check(set(b) == expected_b and not end_differential(pre["delta"], b)
          and not apply_mul(b, b) and not mc_residue(pre["delta"], b)[0],
          "the four-arrow S69 switch is an exact Maurer-Cartan element")
    check(degrees == [Fraction(-8), Fraction(0), Fraction(0), Fraction(8)]
          and all(degree[1] == Fraction(-1) for degree in b_bigradings),
          "the S69 switch is delta-homogeneous but quantum-inhomogeneous")
    check((smooth["long"], smooth["short"], deformed["homology"], smooth["finite"])
          == (5, 4, 9, 9),
          "direct smoothing and all-word algebra both give target dimension 9")

    expected_physical = {
        "S69": {
            (0, ("R", 1), 19), (26, ("R", 1), 27),
            (0, ("R", 1), 27), (26, ("R", 1), 19),
        },
        "S18": {
            (4, ("R", 1), 3), (22, ("R", 1), 23),
            (22, ("R", 1), 3), (4, ("R", 1), 23),
        },
        "S25": {
            (6, ("L", 1), 5), (20, ("L", 1), 21),
            (20, ("L", 1), 5), (6, ("L", 1), 21),
        },
        "S74": {
            (20, ("L", 1), 21), (32, ("L", 1), 33),
            (20, ("L", 1), 33), (32, ("L", 1), 21),
        },
    }
    expected_summands = {
        "S69": [5, 4], "S18": [4, 5],
        "S25": [4, 5], "S74": [7, 2],
    }
    named_orbits = orbit_group(blue)
    physical = {}
    for name, target in targets.items():
        hit = locate_orbit(named_orbits, target)
        if hit is None:
            raise AssertionError(f"missing physical smoothing orbit {name}")
        support = (hit[1], hit[2])
        radius_results = [
            physical_smoothing_pairing(
                red, blue, support, data, pairing=1, eps=eps,
                red_type_d=(red_data, red_vertices, red_delta))
            for eps in (0.003, 0.006, 0.009)
        ]
        result = radius_results[1]
        physical[name] = result
        term_degrees = [
            morphism_term_bigrading(term, grades) for term in result["b"]]
        check(
            result["module_matches"]
            and set(result["b"]) == expected_physical[name]
            and not end_differential(pre["delta"], result["b"])
            and not apply_mul(result["b"], result["b"])
            and all(degree[1] == Fraction(-1) for degree in term_degrees),
            f"{name} smoothing is an explicit four-arrow MC switch")
        check(
            result["wrapped_summands"] == expected_summands[name]
            and result["wrapped"] == result["finite"] == 9,
            f"{name} direct smoothing has exact wrapped rank 9")
        check(all(
            candidate["b"] == result["b"] and candidate["wrapped"] == 9
            for candidate in radius_results),
            f"{name} switch is stable for three surgery radii")

    permutation = {vertex: vertex for vertex in pre["vertices"]}
    permutation.update({4: 22, 22: 4, 5: 21, 21: 5})
    s18_relabelled = collect_f2([
        (permutation[i], word, permutation[j])
        for i, word, j in physical["S18"]["delta"]
    ])
    check(set(s18_relabelled) == set(physical["S25"]["delta"])
          and all(data["unique"][i]["arc"]
                  == data["unique"][permutation[i]]["arc"]
                  for i in pre["vertices"]),
          "S18 and S25 are strictly isomorphic by (4 22)(5 21)")

    def normalized_component_classes(result):
        classes = []
        for index in result["selected"]:
            value = result["classes"][index]
            if value[0] < 0 or (value[0] == 0 and value[1] < 0):
                value = (-value[0], -value[1])
            classes.append(value)
        return tuple(sorted(classes))

    component_classes = {
        name: normalized_component_classes(result)
        for name, result in physical.items()
    }
    check(component_classes == {
        "S69": ((2, 1), (17, 8)),
        "S18": ((2, 2), (17, 6)),
        "S25": ((2, 2), (17, 6)),
        "S74": ((3, 1), (15, 8)),
    }, "the four smoothings have exactly three lifted component-class multisets")
    print("  physical rank-nine summands=" + str({
        name: result["wrapped_summands"] for name, result in physical.items()
    }))
    print("  scope: exact in Tw(B); the instanton tangle-object identification remains open")
    return {
        "base": base["homology"], "deformed": deformed["homology"],
        "b": b, "physical": physical, "component_classes": component_classes,
        "s18_s25_permutation": permutation,
    }


def report_s69_surgery_kwz(red, blue, orbs):
    """Encode the actual S69 smoothing and compare finite and wrapped ranks."""
    from bigons import edges_of, intersections_detailed
    from polygons import self_intersections_detailed
    from surgery_check import (EPS, bigon_matrix_multi, iota_pairing,
                               smooth_curve)

    print("== exact KWZ pairing of the S69 geometric smoothing ==")
    hit = locate_orbit(orbs, TARGETS["S69"])
    if hit is None:
        check(False, "S69 is available for the surgery certificate")
        return None
    _, _, preimages = hit
    c1, c2 = preimages
    p2, mismatch = iota_pairing(edges_of(blue), c1, 1, c2, EPS)
    comps = smooth_curve(blue, [(c1, 1), (c2, p2)], EPS)
    classes = [deck_of_path(comp) for comp in comps]
    print(f"  torus components={len(comps)}; classes={classes}")
    check(p2 == 1 and mismatch < 1e-8,
          "pairing 1 is iota-compatible at the two S69 lifts")
    check(classes == [(17, 8), (-2, -1), (-2, -1)],
          "S69 smoothing has torus classes (17,8),(-2,-1),(-2,-1)")

    red_data, red_vertices, red_delta = _q7_red_type_d(red)
    encoded = []
    for comp in comps:
        comp_data = encode_type_d(comp)
        if len(comp_data["hits"]) == len(comp_data["unique"]):
            comp_pre = cyclic_precurve(comp_data)
        else:
            comp_pre = geometric_precurve(comp_data)
        encoded.append((comp_data, comp_pre))

    long_index = classes.index((17, 8))
    short_indices = [i for i, value in enumerate(classes) if value == (-2, -1)]
    long_data, long_pre = encoded[long_index]
    matching_shorts = [
        i for i in short_indices
        if encoded[i][1]["vertices"] == red_vertices
        and encoded[i][1]["delta"] == red_delta
    ]
    check(len(long_data["hits"]) == 60
          and len(long_data["unique"]) == 30
          and len(long_pre["vertices"]) == 24
          and len(long_pre["cancellations"]) == 3
          and not long_pre["residue"],
          "the long projected component cancels to a 24-dot type-D object")
    check(len(matching_shorts) == 1,
          "one oriented short lift has exactly the red earring type-D cycle")

    long_self = len(self_intersections_detailed(comps[long_index]))
    long_short = [
        len(intersections_detailed(comps[long_index], comps[i]))
        for i in short_indices
    ]
    print(f"  residual torus crossings: long self={long_self}; "
          f"long/short={long_short}")
    check(long_self == 102 and long_short == [29, 29],
          "the smoothing is not in the embedded/no-intersection surgery regime")

    finite_gens, finite_delta = bigon_matrix_multi(red, comps)
    finite_rank = rank_f2(finite_delta)
    finite_homology = len(finite_gens) - 2 * finite_rank
    check(len(finite_gens) == 13 and finite_rank == 2
          and finite_homology == 9,
          "the finite single-traversal smoothing statistic is 9")

    long_result = exact_red_blue_homology(
        red_data, red_vertices, red_delta,
        long_data, long_pre, long_pre["delta"])
    red_end = exact_red_blue_homology(
        red_data, red_vertices, red_delta,
        red_data, {"vertices": red_vertices}, red_delta)
    wrapped_total = long_result["homology"] + red_end["homology"]
    print(f"  exact wrapped dimensions: red-long={long_result['homology']}; "
          f"red-red={red_end['homology']}; total={wrapped_total}")
    check(long_result["homology"] == 4
          and red_end["homology"] == 4
          and wrapped_total == 8,
          "the exact wrapped pairing of the S69 smoothing is 4+4=8, not 9")
    return {
        "classes": classes,
        "finite": finite_homology,
        "long": long_result["homology"],
        "short": red_end["homology"],
        "wrapped": wrapped_total,
    }


def _square_zero_masks(representatives):
    """All masks whose represented sum has square zero, including zero."""
    n = len(representatives)
    products = {}
    raw = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            raw[i][j] = apply_mul(representatives[i], representatives[j])
            for term in raw[i][j]:
                products.setdefault(term, len(products))

    def product_bits(terms):
        value = 0
        for term in terms:
            value ^= 1 << products[term]
        return value

    diagonal = [product_bits(raw[i][i]) for i in range(n)]
    cross = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            cross[i][j] = cross[j][i] = product_bits(
                collect_f2(raw[i][j] + raw[j][i]))
    solutions = [0]
    previous_gray = 0
    square = 0
    interaction = [0] * n
    for integer in range(1, 1 << n):
        gray = integer ^ (integer >> 1)
        changed = gray ^ previous_gray
        i = changed.bit_length() - 1
        square ^= diagonal[i] ^ interaction[i]
        for j in range(n):
            if j != i:
                interaction[j] ^= cross[i][j]
        if square == 0:
            solutions.append(gray)
        previous_gray = gray
    return solutions


def _face_tail_homology(face, red_data, red_vertices, red_delta,
                         blue_data, blue_pre, blue_delta):
    """Exact F2-dimension of one positive face tail in the U-linear cases."""
    labels, matrix = _face_tower_complex(
        face, red_data, red_vertices, red_delta,
        blue_data, blue_pre, blue_delta)
    active, residual, _, _ = _cancel_unit_tower_pairs(
        labels, dict(matrix), [])
    if any(coeff != 2 for coeff in residual.values()):
        raise AssertionError("face tail is not U times a constant matrix")
    order = {old: new for new, old in enumerate(sorted(active))}
    columns = [[] for _ in active]
    for (row, col), coeff in residual.items():
        if coeff == 2:
            columns[order[col]].append(order[row])
    rank = sparse_column_rank(columns)
    if 2 * rank != len(active):
        raise AssertionError("face tail has free F2[U] homology")
    for column in columns:
        image = 0
        for middle in column:
            for row in columns[middle]:
                image ^= 1 << row
        if image:
            raise AssertionError("face-tail constant matrix does not square to zero")
    return rank


def report_relative_pairing_probe(red, data):
    """Full red-blue Hom homology, with finite-support cross-checks."""
    print("== exact red-blue Hom calculation over F2[U] ==")
    red_data, red_vertices, red_delta = _q7_red_type_d(red)
    pre = data["precurve"]
    u0 = [(7, ("M", 2), 12), (4, ("M", 2), 15)]
    ustar = [(1, ("M", 1), 27), (0, ("M", 2), 24)]
    uaux = [(1, ("M", 1), 33)]
    ucombo = collect_f2(ustar + uaux)
    deformations = {
        "base": [], "u0": u0, "a": uaux,
        "u*": ustar, "u*+a": ucombo,
    }
    exact = {
        label: exact_red_blue_homology(
            red_data, red_vertices, red_delta, data, pre,
            collect_f2(pre["delta"] + cycle))
        for label, cycle in deformations.items()
    }
    expected_faces = {
        "L": {"tower_rank": 20, "unit_pairs": 0,
              "residual_rank": 20, "homology": 10},
        "M": {"tower_rank": 180, "unit_pairs": 40,
              "residual_rank": 100, "homology": 50},
        "R": {"tower_rank": 80, "unit_pairs": 0,
              "residual_rank": 80, "homology": 40},
    }
    check(all(result["identities"] == 100
              and result["tail_homology"] == 100
              and result["faces"] == expected_faces
              for result in exact.values()),
          "positive towers reduce to 100 copies of F2[U]/(U)")
    exact_dims = {label: result["homology"]
                  for label, result in exact.items()}
    beta_ranks = {label: result["beta_rank"]
                  for label, result in exact.items()}
    print(f"  identity generators=100; positive-tail H=10+50+40=100")
    print(f"  identity-to-tail ranks: {beta_ranks}")
    print(f"  full Hom dimensions: {exact_dims}")
    check(beta_ranks == {
        "base": 97, "u0": 98, "a": 97, "u*": 96, "u*+a": 96,
    } and exact_dims == {
        "base": 6, "u0": 4, "a": 6, "u*": 8, "u*+a": 8,
    }, "mapping-cone formula gives exact dimensions 6,4,6,8,8")

    print("  -- independent bounded-support cross-check --")
    dimensions = {label: [] for label in deformations}
    for cutoff in (1, 2, 4, 8):
        row = {}
        for label, cycle in deformations.items():
            blue_delta = collect_f2(pre["delta"] + cycle)
            samples = [finite_support_homology(
                red_data, red_vertices, red_delta, data, pre, blue_delta,
                cutoff, primitive_extra=extra)
                for extra in (0, 2, 4)]
            check(samples[0] == samples[1] == samples[2],
                  f"{label} boundary space stabilizes at cutoff {cutoff}")
            row[label] = samples[0]
            dimensions[label].append(samples[0]["homology"])
        print(f"  cutoff {cutoff}: chain={row['base']['chain']}; "
              + ", ".join(
                  f"{label} H={row[label]['homology']}"
                  for label in deformations))
    check(dimensions == {
        "base": [6, 6, 6, 6], "u0": [4, 4, 4, 4],
        "a": [6, 6, 6, 6], "u*": [8, 8, 8, 8],
        "u*+a": [8, 8, 8, 8],
    }, "bounded-support dimensions agree with the exact polynomial result")

    base = exact_dims["base"]
    plus_two_cycles = []
    rank_nine_cycles = []
    tested_mc = 0
    tested_non_mc = 0
    for qdegree, reps in end_h_representatives().items():
        for mask in range(1, 1 << len(reps)):
            cycle = collect_f2([
                term for k, rep in enumerate(reps) if mask & (1 << k)
                for term in rep])
            if apply_mul(cycle, cycle):
                tested_non_mc += 1
                continue
            tested_mc += 1
            result = exact_red_blue_homology(
                red_data, red_vertices, red_delta, data, pre,
                collect_f2(pre["delta"] + cycle))
            if result["homology"] - base == 2:
                plus_two_cycles.append((qdegree, cycle))
            if result["homology"] == 9:
                rank_nine_cycles.append((qdegree, cycle))
    normalize = lambda cycle: tuple(sorted(cycle, key=repr))
    plus_two_set = {normalize(cycle) for _, cycle in plus_two_cycles}
    check(len(plus_two_cycles) == 2
          and {qdegree for qdegree, _ in plus_two_cycles} == {Fraction(-8)}
          and plus_two_set == {normalize(ustar), normalize(ucombo)},
          "the homogeneous MC screen isolates the +2 coset {u*,u*+a}")
    check(not rank_nine_cycles,
          "no homogeneous Maurer-Cartan class gives Hom dimension 9")
    print(f"  homogeneous cohomology screen: {tested_mc} MC classes, "
          f"{tested_non_mc} classes with nonzero square; "
          "+2 shift exactly {u*,u*+a}; no tested class has dimension 9")

    # Forgetting quantum degree permits sums across the homogeneous sectors.
    # All displayed representatives live in a single face, and products
    # between distinct faces vanish, so the square-zero equation splits into
    # independent L and M equations.
    all_reps = [rep for qdegree in sorted(end_h_representatives())
                for rep in end_h_representatives()[qdegree]]
    face_reps = {"L": [], "M": []}
    for rep in all_reps:
        faces = {word[0] for _, word, _ in rep}
        if len(faces) != 1 or next(iter(faces)) not in face_reps:
            raise AssertionError("End representative does not lie in one tested face")
        face_reps[next(iter(faces))].append(rep)
    face_masks = {
        face: _square_zero_masks(reps) for face, reps in face_reps.items()
    }
    check(len(face_masks["L"]) == 4 and len(face_masks["M"]) == 18432,
          "ungraded MC equation splits into 4 L-face and 18432 M-face solutions")
    face_tail_values = {}
    for face, masks in face_masks.items():
        values = set()
        for mask in masks:
            cycle = collect_f2([
                term for k, rep in enumerate(face_reps[face])
                if mask & (1 << k) for term in rep
            ])
            values.add(_face_tail_homology(
                face, red_data, red_vertices, red_delta, data, pre,
                collect_f2(pre["delta"] + cycle)))
        face_tail_values[face] = values
    check(face_tail_values == {"L": {10}, "M": {50}},
          "all 18436 facewise MC cases preserve tail homology 10+50+40=100")
    ungraded_mc = len(face_masks["L"]) * len(face_masks["M"]) - 1
    check(ungraded_mc == 73727,
          "the displayed 22-class span contains 73727 nonzero ungraded MC sums")
    check((100 + 100) % 2 == 0,
          "every such ungraded deformation has even full Hom dimension, never 9")
    print("  exhaustive ungraded span: 73,727 nonzero MC sums; all have "
          "100-dimensional positive tail, hence even full Hom dimension")
    print("  scope warning: the algebraic all-word calculation is exact for the "
          "encoded type-D objects; the character-variety/tangle/instanton "
          "identification remains open")
    return exact_dims


def report_s69_resolution(blue):
    print("== KWZ resolution of S69 (Def. 5.17, right-turn paths) ==")
    rec = named_s69_chart(blue)
    if rec is None or rec["xy"] is None or rec["face"] != "R":
        check(False, "S69 sits in the right face of the chart")
        return None
    _, orbs = orbit_records(blue)
    hit = locate_orbit(orbs, TARGETS["S69"])
    _, _, pre = hit
    res = resolve_crossing_in_face(
        blue, pre, "R", rec["xy"], PUNCTURE["R"], xc=0.5)
    print(f"  carriers found: {res['carriers']}")
    print(f"  kept right-turn terms: {len(res['terms'])}")
    for tag, wind, a, b in res["terms"]:
        print(f"    {tag}  D^{wind}  {tuple(round(v, 3) for v in a)} -> "
              f"{tuple(round(v, 3) for v in b)}")
    check(res["carriers"] == 2, "exactly two charted branches carry S69")
    crosses = [t for t in res["terms"]
               if math.hypot(t[2][0] - t[3][0], t[2][1] - t[3][1]) > 0.05]
    print(f"  distinct-endpoint terms (candidate φ(S69) in A_R): {len(crosses)}")
    return res


# ---------------------------------------------------------------------------
# reports
# ---------------------------------------------------------------------------
def report_records(blue, named):
    print("== q=7 physical records at default parameters ==")
    print(f"blue edges: {len(blue) - 1}")
    for name in ("S69", "S18", "S25", "S74"):
        rec = named[name]
        if rec is None:
            print(f"  {name}: NOT FOUND")
            continue
        print(f"  {name} = code {rec['name']}  P={tuple(round(v, 6) for v in rec['P'])}  "
              f"dist={rec['dist']:.3e}  preimages={rec['n_pre']}")
        for p in rec["preimages"]:
            print(f"    pre{p['index']}: edges ({p['kA']},{p['kB']})  "
                  f"t=({p['tA']:.6f},{p['tB']:.6f})  "
                  f"deg A->B={p['degrees']['A->B']} B->A={p['degrees']['B->A']}  "
                  f"deck A->B={p['deck_A_to_B']} B->A={p['deck_B_to_A']}  "
                  f"short={p['deck_short']}  "
                  f"param=({p['paramA']:.6f},{p['paramB']:.6f})")


def report_witness(blue, orbs, named):
    print("== S69 straight-through self-bigon witness ==")
    rec = named["S69"]
    if rec is None or rec["n_pre"] < 2:
        check(False, "S69 has two preimages")
        return None
    # Handoff: preimage 1 has edges (845, 941).
    pre = rec["preimages"]
    pre1 = min(pre, key=lambda p: abs(p["kB"] - 941) + abs(p["kA"] - 845))
    print(f"  S69 preimage-1 edges ({pre1['kA']},{pre1['kB']})  "
          f"paramB={pre1['paramB']:.9f}")
    w = find_handoff_witness(blue, orbs, pre1)
    if w is None:
        check(False, "found crossings with handoff edge pairs (6,905) and (43,945)")
        return None
    print(f"  corner0 S{w['idx0']} edges ({w['s0']['kA']},{w['s0']['kB']})")
    print(f"  corner1 S{w['idx1']} edges ({w['s1']['kA']},{w['s1']['kB']})")
    if w["typed"] is None:
        check(False, "accepted self-bigon on those two crossings")
        return w
    print(f"  typing sel={w['typed']['sel']}  orients={w['typed']['fbs']}")
    print(f"  legs through S69-pre1 param {pre1['paramB']:.9f}: {w['hits']}")
    check(bool(w["hits"]),
          "one bigon arc passes through the S69 preimage-1 point "
          "without a branch change")
    return w


def report_chart():
    print("== pillowcase-to-disk chart at (0,0) ==")
    imgs = chart_corner_images()
    print("  corner images:", {k: (None if v is None else (round(v[0], 6), round(v[1], 6)))
                               for k, v in imgs.items()})
    check(imgs[(0.0, 0.0)] is None, "(0,0) is the projection centre")
    others = [imgs[(PI, 0.0)], imgs[(0.0, PI)], imgs[(PI, PI)]]
    check(all(p is not None for p in others), "the other three corners are finite")
    if all(p is not None for p in others):
        xs = sorted(p[0] for p in others)
        ys = [abs(p[1]) for p in others]
        print(f"  other x-coords (sorted) = {[round(x, 6) for x in xs]}")
        print(f"  other |y| = {[round(y, 8) for y in ys]}")
        check(all(y < 1e-9 for y in ys), "other corners land on the x-axis")
        # Handoff: they map to (-1,0), (0,0), (1,0).
        check(abs(xs[0] + 1) < 1e-6 and abs(xs[1]) < 1e-6 and abs(xs[2] - 1) < 1e-6,
              "other corners map to (-1,0), (0,0), (1,0)")


def report_finite_complex(red, blue, named):
    print("== default finite bigon matrix ==")
    gens, d = bigon_matrix(red, blue)
    ents = [(i, j) for i in range(len(d)) for j in range(len(d)) if d[i][j]]
    rk = rank_f2(d)
    print(f"  gens={len(gens)}  entries={ents}  rank={rk}  h={len(gens) - 2 * rk}")
    check(len(gens) == 13, "13 red-blue generators")
    check(ents == [(0, 11), (6, 4), (7, 5)], "bigon entries (0,11),(6,4),(7,5)")
    check(rk == 3 and len(gens) - 2 * rk == 7, "rank 3 and statistic 7")
    rec = named["S69"]
    if rec is not None:
        check(rec["n_pre"] == 2, "S69 has two torus preimages")
        edges = sorted((p["kA"], p["kB"]) for p in rec["preimages"])
        print(f"  S69 edge pairs: {edges}")
        check(edges == [(522, 618), (845, 941)] or
              sorted(tuple(sorted(e)) for e in edges) ==
              sorted(tuple(sorted(e)) for e in [(522, 618), (845, 941)]),
              "S69 edge pairs are (522,618) and (845,941)")
        pre0 = next(p for p in rec["preimages"] if (p["kA"], p["kB"]) == (522, 618))
        pre1 = next(p for p in rec["preimages"] if (p["kA"], p["kB"]) == (845, 941))
        check(pre0["degrees"]["B->A"] == 1 and pre1["degrees"]["A->B"] == 1,
              "S69 degree-one jumps are pre0 B->A and pre1 A->B")
        check(pre0["deck_short"] == (-2, -1) and pre1["deck_short"] == (-2, -1),
              "S69 short deck class is ±(2,1)")
        check(abs(pre1["paramB"] - WITNESS_THROUGH_PARAM) < 1e-5,
              "S69 preimage-1 parameter is 941.640214...")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    only_w = "--witness" in argv
    only_r = "--records" in argv
    only_e = "--encode" in argv
    print("== q7_kwz: records, witness, chart, type-D (Q7_HANDOFF_2026-08-12) ==\n")
    red, blue, _ = build_q7()
    named, recs, orbs = named_records(blue)
    if not only_w and not only_e:
        report_finite_complex(red, blue, named)
        print()
        report_records(blue, named)
        print()
        print(f"== KWZ special puncture: original pillowcase corner {KWZ_SPECIAL} ==")
        check(min(_tdist(point, KWZ_SPECIAL) for point in blue[:-1]) > 0.1,
              "the q=7 blue arc misses the chosen special puncture")
        report_chart()
        print()
    if not only_r and not only_e:
        report_witness(blue, orbs, named)
        print()
    if not only_w and not only_r:
        red_kwz = translate_to_kwz_chart(red)
        blue_kwz = translate_to_kwz_chart(blue)
        targets_kwz = {
            name: translate_target_to_kwz(point)
            for name, point in TARGETS.items()
        }
        data = report_encoding(blue_kwz, s69_target=targets_kwz["S69"])
        print()
        report_admissible_kwz(red_kwz, blue_kwz, targets_kwz, data)
        print()
    if FAILS:
        print("FAILURES: " + "; ".join(FAILS))
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
