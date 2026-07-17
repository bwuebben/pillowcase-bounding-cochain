#!/usr/bin/env python3
"""
resolve.py -- step (b): the perturbed tangle sum (RESEARCH_LOG sec 25).

Implements Smith's Theorem thm:cross (arXiv:2412.06066, read in TeX source
2026-07-16), whose content we independently pinned in secs 23-24: perturbing the
middle piece C_3 of a Conway sum resolves each seam fiber circle by CUT-AND-PASTE:

  * Each fiber circle (seam s, arc [lo,hi]) meets the abelian locus at its two fold
    endpoints rho_0 (theta=lo) and rho_pi (theta=hi) -- cone-on-4-points singular
    points of the unperturbed variety.
  * The resolution deletes the circle and reconnects the four main-curve branch
    ends DIAGONALLY:  A+(lo) <-> A-(hi)  and  A-(lo) <-> A+(hi),  where A+ is the
    branch on the (0,pi)-annulus side of the seam and A- the (pi,2pi) side.
    (The s-sign data of Smith's lem:endop only routes the invisible internal
    halves; the visible outcome is always this diagonal.)
  * Each connector crosses the seam once; the two connectors of one circle cross
    each other once => ONE new P self-intersection per resolved circle (support
    candidates for the bounding cochain), plus transversal crossings with whatever
    other strands come near the seam arc.

Validated by hand + figure against Smith's worked example R_t(Q_{1/2}+Q_{-1/3})
(his Figures Ex3/ExResA/ExResB): one circle at gamma=pi with arc [pi/6, 5pi/6];
the resolution splits the curve into TWO components -- a corner-to-corner arc
through 2 fundamental-domain strands and a closed loop through 4.

Everything lives on the double cover T^2 = [0,2pi)^2 with iota(g,t)=(-g,-t);
P self-crossings = T^2 self-crossings / 2 (iota-orbits).
"""
import math
from tangles import (TAU, PI, curve, segments, tangle_sum, fiber_circles,
                     south_twists, east_twists)

TOL = 1e-6


def _norm(p):
    return (p[0] % TAU, p[1] % TAU)


def _tdist(p, q):
    dg = abs(p[0] % TAU - q[0] % TAU)
    dt = abs(p[1] % TAU - q[1] % TAU)
    return math.hypot(min(dg, TAU - dg), min(dt, TAU - dt))


class _EndIndex:
    """Spatial hash of segment endpoints on T^2 with tolerance lookup."""

    def __init__(self):
        self.buckets = {}

    def _bk(self, p):
        return (int((p[0] % TAU) / 1e-4), int((p[1] % TAU) / 1e-4))

    def add(self, p, item):
        self.buckets.setdefault(self._bk(p), []).append((p, item))

    def pop_near(self, p, tol=TOL):
        bx, by = self._bk(p)
        nb = int(TAU / 1e-4)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                key = ((bx + dx) % (nb + 1), (by + dy) % (nb + 1))
                lst = self.buckets.get(key, [])
                for k, (q, item) in enumerate(lst):
                    if _tdist(p, q) < tol:
                        lst.pop(k)
                        return item
        return None


def assemble(segs):
    """Chain a segment soup (g0,g1,t0,t1) into polylines on T^2 by tolerance-based
    endpoint matching. Returns list of polylines (lists of (gamma,theta))."""
    segs = list(segs)
    idx = _EndIndex()
    for i, (g0, g1, t0, t1) in enumerate(segs):
        idx.add((g0, t0), (i, 0))
        idx.add((g1, t1), (i, 1))
    used = [False] * len(segs)

    def take(i):
        used[i] = True
        g0, g1, t0, t1 = segs[i]
        # remove both endpoint records of segment i from the index
        for p, e in (((g0, t0), 0), ((g1, t1), 1)):
            bx = idx._bk(p)
            lst = idx.buckets.get(bx, [])
            for k, (q, item) in enumerate(lst):
                if item == (i, e):
                    lst.pop(k)
                    break

    polylines = []
    for i0 in range(len(segs)):
        if used[i0]:
            continue
        take(i0)
        g0, g1, t0, t1 = segs[i0]
        pts = [(g0, t0), (g1, t1)]
        # extend forward from the tail, then backward from the head
        for forward in (True, False):
            while True:
                tip = pts[-1] if forward else pts[0]
                hit = idx.pop_near(tip)
                if hit is None:
                    break
                j, e = hit
                if used[j]:
                    continue
                take(j)
                a, b = (segs[j][0], segs[j][2]), (segs[j][1], segs[j][3])
                nxt = b if e == 0 else a
                if forward:
                    pts.append(nxt)
                else:
                    pts.insert(0, nxt)
        polylines.append(pts)
    return polylines


def is_closed(poly):
    return _tdist(poly[0], poly[-1]) < TOL


# ---------------------------------------------------------------------------
# the cut-and-paste resolution
# ---------------------------------------------------------------------------
def _connector_U(s, lo, hi, side_a, eps, zero_at, N=160):
    """Connector from (s,lo) to (s,hi), leaving lo on side_a (+1 = the (0,pi)
    annulus, -1 = the (pi,2pi) annulus), arriving at hi on -side_a; crosses the
    seam once, at fraction zero_at."""
    sgn_I = 1.0 if abs(s % TAU) < 1e-9 else -1.0   # gamma-offset sign giving side (0,pi)
    pts = []
    for k in range(N + 1):
        w = k / N
        if w <= zero_at:
            prof = math.sin(PI * w / zero_at)
        else:
            prof = -math.sin(PI * (w - zero_at) / (1 - zero_at))
        off = side_a * sgn_I * eps * prof
        pts.append(((s + off) % TAU, (lo + (hi - lo) * w) % TAU))
    return pts


def _iota(poly):
    return [((-g) % TAU, (-t) % TAU) for (g, t) in poly]


def _side_of_point(p, s):
    """Side of T^2 point p relative to seam s: +1 if gamma in (0,pi), else -1."""
    g = p[0] % TAU
    return 1.0 if 0 < g < PI else -1.0


def _cut(polylines, cutpts):
    """Split polylines at the given T^2 seam points. Returns open chains; every
    chain end at a cut point is tagged. Assumes cut points are existing vertices."""
    chains = []
    for poly in polylines:
        hit = []
        for vi, p in enumerate(poly[:-1] if is_closed(poly) else poly):
            for (s, th) in cutpts:
                if _tdist(p, (s, th)) < 1e-5:
                    hit.append(vi)
                    break
        hit = sorted(set(hit))
        if not hit:
            chains.append(poly)
            continue
        if is_closed(poly):
            core = poly[:-1]
            n = len(core)
            base = hit[0]
            rot = core[base:] + core[:base] + [core[base]]
            offs = sorted([(h - base) % n for h in hit]) + [n]
            for k in range(len(offs) - 1):
                piece = rot[offs[k]:offs[k + 1] + 1]
                if len(piece) > 1:
                    chains.append(piece)
        else:
            offs = sorted(set([0] + hit + [len(poly) - 1]))
            for k in range(len(offs) - 1):
                piece = poly[offs[k]:offs[k + 1] + 1]
                if len(piece) > 1:
                    chains.append(piece)
    return chains


def resolve(segs, circles, eps=0.05):
    """Perturbed tangle sum (thm:cross cut-and-paste). Returns (polylines, xinfo):
    the resolved immersed curve on T^2, and per-strip connector records."""
    polylines = assemble(segs)
    # per circle: strip-U pair (lo,hi) and strip-D pair (2pi-hi, 2pi-lo)
    strips = []
    for ci, c in enumerate(circles):
        s, lo, hi = c['seam'], c['lo'], c['hi']
        strips.append((s, lo, hi, ci, 'U'))
        strips.append((s, (TAU - hi) % TAU, (TAU - lo) % TAU, ci, 'D'))
    cutpts = []
    for (s, a, b, ci, tag) in strips:
        cutpts += [(s, a), (s, b)]
    chains = _cut(polylines, cutpts)

    # port registry: (cut point, side) -> chain end
    ports = {}
    for chi, ch in enumerate(chains):
        for endpos, inward in ((0, 1), (len(ch) - 1, -2)):
            p = ch[endpos]
            for cp in cutpts:
                if _tdist(p, cp) < 1e-5:
                    side = _side_of_point(ch[endpos + (1 if endpos == 0 else -1)], cp[0])
                    ports.setdefault((round(cp[0], 6), round(cp[1] % TAU, 6), side),
                                     []).append((chi, endpos))
                    break

    # build connectors: strip U explicitly, strip D as exact iota-image of U
    connectors, xinfo = [], []
    for ci, c in enumerate(circles):
        s, lo, hi = c['seam'], c['lo'], c['hi']
        eps_c = eps * (1.0 + 0.35 * (ci % 3))
        u1 = _connector_U(s, lo, hi, +1.0, eps_c, zero_at=0.42)
        u2 = _connector_U(s, lo, hi, -1.0, eps_c * 0.62, zero_at=0.58)
        d1, d2 = _iota(u1), _iota(u2)
        # ends: u1: (lo,+1)->(hi,-1); u2: (lo,-1)->(hi,+1)
        # d1 = iota(u1): (2pi-lo,-1)->(2pi-hi,+1); d2: (2pi-lo,+1)->(2pi-hi,-1)
        la, lb = (TAU - lo) % TAU, (TAU - hi) % TAU
        connectors.append((u1, (s, lo, +1.0), (s, hi, -1.0)))
        connectors.append((u2, (s, lo, -1.0), (s, hi, +1.0)))
        connectors.append((d1, (s, la, -1.0), (s, lb, +1.0)))
        connectors.append((d2, (s, la, +1.0), (s, lb, -1.0)))
        xinfo.append(dict(circle=ci, seam=s, U=(u1, u2), D=(d1, d2)))

    # stitch: walk chains and connectors through the ports
    def port_key(pt, side):
        return (round(pt[0], 6), round(pt[1] % TAU, 6), side)

    # map each (point,side) to its unique chain end and unique connector end
    con_ends = {}
    for coni, (cpts, endA, endB) in enumerate(connectors):
        for (pt, which) in ((endA, 0), (endB, 1)):
            key = port_key((pt[0], pt[1]), pt[2])
            assert key not in con_ends, f"duplicate connector port {key}"
            con_ends[key] = (coni, which)

    chain_ends = {}
    for key, lst in ports.items():
        assert len(lst) == 1, f"port {key} has {len(lst)} chain ends (expected 1)"
        chain_ends[key] = lst[0]

    assert set(con_ends.keys()) == set(chain_ends.keys()), (
        "connector ports and chain-end ports mismatch:\n"
        f"  only-connector: {sorted(set(con_ends) - set(chain_ends))[:4]}\n"
        f"  only-chain:     {sorted(set(chain_ends) - set(con_ends))[:4]}")

    used_ch, used_con = set(), set()
    out = []
    for chi0 in range(len(chains)):
        if chi0 in used_ch:
            continue
        # closed untouched chains pass through
        ch = chains[chi0]
        end_key = None
        for key, (c_i, pos) in chain_ends.items():
            if c_i == chi0:
                end_key = key
                break
        if end_key is None:
            used_ch.add(chi0)
            out.append(ch)
            continue
        # walk: chain -> connector -> chain -> ... until loop closes
        path = []
        cur_chain = chi0
        # orient so we END at a port (start from the other end)
        (c_i, pos) = chain_ends[end_key]
        pts = list(chains[cur_chain])
        if pos == 0:
            pts = pts[::-1]  # make the ported end the tail
        used_ch.add(cur_chain)
        path = pts
        while True:
            tail = path[-1]
            # find the port at the tail
            side = None
            for key, (c_i, pos) in chain_ends.items():
                if c_i == cur_chain and _tdist(chains[cur_chain][pos], tail) < 1e-9:
                    pass
            # identify port by location+side of the path's last inward vertex
            s_here = None
            for cp in cutpts:
                if _tdist(tail, cp) < 1e-5:
                    s_here = cp
                    break
            if s_here is None:
                break  # open end not at a port (e.g. numerical corner) -- stop
            side = _side_of_point(path[-2], s_here[0])
            key = port_key(s_here, side)
            if key not in con_ends:
                break
            coni, which = con_ends[key]
            if coni in used_con:
                break  # closed the loop
            used_con.add(coni)
            cpts_, endA, endB = connectors[coni]
            cpts_ = list(cpts_) if which == 0 else list(cpts_)[::-1]
            path += cpts_[1:]
            # continue with the chain plugged into the connector's far end
            far = (endB if which == 0 else endA)
            fkey = port_key((far[0], far[1]), far[2])
            nchi, npos = chain_ends[fkey]
            if nchi in used_ch:
                # arrived back at the start chain: close up
                break
            used_ch.add(nchi)
            npts = list(chains[nchi])
            if npos != 0:
                npts = npts[::-1]
            path += npts[1:]
            cur_chain = nchi
        out.append(path)
    return out, xinfo


# ---------------------------------------------------------------------------
# torus-aware intersection counting
# ---------------------------------------------------------------------------
def _seg_int_torus(p1, p2, p3, p4):
    mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
    bx, by = (p3[0] + p4[0]) / 2, (p3[1] + p4[1]) / 2
    sx = TAU * round((mx - bx) / TAU)
    sy = TAU * round((my - by) / TAU)
    x1, y1 = p1; x2, y2 = p2
    x3, y3 = p3[0] + sx, p3[1] + sy
    x4, y4 = p4[0] + sx, p4[1] + sy
    d = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3)
    if abs(d) < 1e-14:
        return None
    t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / d
    u = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / d
    if 1e-9 < t < 1 - 1e-9 and 1e-9 < u < 1 - 1e-9:
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))
    return None


def _short_edges(poly):
    """Edges with locally-unwrapped second endpoint (edges are short)."""
    out = []
    for k, (p, q) in enumerate(zip(poly, poly[1:])):
        g0, t0 = p[0] % TAU, p[1] % TAU
        dg = (q[0] - p[0] + PI) % TAU - PI
        dt = (q[1] - p[1] + PI) % TAU - PI
        out.append((k, (g0, t0), (g0 + dg, t0 + dt)))
    return out


def self_crossings_T2(polylines):
    """Transverse self-intersections of the curve on T^2 (adjacent edges of the
    same polyline excluded)."""
    E = []
    for pi_, poly in enumerate(polylines):
        n = len(poly) - 1
        closed = is_closed(poly)
        for (k, a, b) in _short_edges(poly):
            E.append((pi_, k, n, closed, a, b))
    pts = []
    for i in range(len(E)):
        pi1, k1, n1, cl1, a1, b1 = E[i]
        for j in range(i + 1, len(E)):
            pi2, k2, n2, cl2, a2, b2 = E[j]
            if pi1 == pi2:
                dk = abs(k1 - k2)
                if dk <= 1 or (cl1 and dk == n1 - 1):
                    continue
            r = _seg_int_torus(a1, b1, a2, b2)
            if r:
                pts.append((r[0] % TAU, r[1] % TAU))
    return pts


def fold_crossings(pts, tol=5e-3):
    """Group T^2 crossings into iota-orbits; P self-crossings = number of orbits."""
    orbits = []
    for p in pts:
        placed = False
        for orb in orbits:
            for q in orb:
                ip = ((-q[0]) % TAU, (-q[1]) % TAU)
                if _tdist(p, q) < tol or _tdist(p, ip) < tol:
                    orb.append(p)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            orbits.append([p])
    return orbits


def P_components(polylines):
    """Number of P-components (iota-invariant T^2 components count once,
    iota-paired ones once per pair)."""
    def invariant(poly):
        sample = poly[:: max(1, len(poly) // 16)]
        for p in sample:
            ip = ((-p[0]) % TAU, (-p[1]) % TAU)
            if not any(_tdist(ip, q) < 8e-2 for q in poly):
                return False
        return True
    inv = sum(1 for p in polylines if invariant(p))
    return inv + (len(polylines) - inv) // 2


# ---------------------------------------------------------------------------
# battery (sec 25)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    results = []

    def check(name, ok):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    print("== resolution of the perturbed tangle sum (RESEARCH_LOG sec 25) ==\n")

    print("(12) Smith's worked example: R_t(Q(1/2)+Q(-1/3))  [his Figs Ex3/ExResA/ExResB]")
    sA = segments(curve(south_twists(2)))[0]
    sB = segments(curve(south_twists(-3)))[0]
    ssum = tangle_sum(sA, sB)
    circ = fiber_circles(sA, sB)
    check(f"one fiber circle at gamma=pi (got {len(circ)})",
          len(circ) == 1 and abs(circ[0]['seam'] - PI) < 1e-9)
    check(f"arc = [pi/6, 5pi/6] (got [{circ[0]['lo']:.4f},{circ[0]['hi']:.4f}])",
          abs(circ[0]['lo'] - PI / 6) < 1e-6 and abs(circ[0]['hi'] - 5 * PI / 6) < 1e-6)
    unres = assemble(ssum)
    check(f"unperturbed sum assembles into 1 closed T^2 component (got {len(unres)})",
          len(unres) == 1 and is_closed(unres[0]))
    res, xd = resolve(ssum, circ)
    ncomp = P_components(res)
    check(f"resolved curve has 2 P-components (got {ncomp})  [= ExResA + ExResB]",
          ncomp == 2)
    xt2 = self_crossings_T2(res)
    orbs = fold_crossings(xt2)
    print(f"        T^2 self-crossings: {len(xt2)}; P self-crossings: {len(orbs)}")
    check("every T^2 crossing has an iota-partner (orbits of size 2)",
          len(xt2) % 2 == 0 and all(len(o) == 2 for o in orbs))
    check("P self-crossings = 3 (1 connector X + 2 connector-vs-R2-bounce)",
          len(orbs) == 3)

    print("\n(13) the BLUE curve of P(-2,3,5): R_t(Q(1/3)+Q(1/5))  [Smith Fig 2bigons]")
    s3 = segments(curve(south_twists(3)))[0]
    s5 = segments(curve(south_twists(5)))[0]
    bsum = tangle_sum(s3, s5)
    bcirc = fiber_circles(s3, s5)
    check(f"4 fiber circles, 2 per seam (got {len(bcirc)})", len(bcirc) == 4)
    check("unperturbed blue assembles into 1 closed T^2 component",
          len(assemble(bsum)) == 1)
    bres, bxd = resolve(bsum, bcirc)
    bxt2 = self_crossings_T2(bres)
    borbs = fold_crossings(bxt2)
    check("iota-pairing of blue crossings clean (orbits of size 2)",
          len(bxt2) % 2 == 0 and all(len(o) == 2 for o in borbs))
    ncomp_b = P_components(bres)
    print(f"        blue: {len(bres)} T^2 components, {ncomp_b} P-components; "
          f"{len(borbs)} P self-crossings ({len(bxt2)} on T^2)")
    nX = sum(1 for x in bxd
             if len(self_crossings_T2([x['U'][0], x['U'][1]])) == 1
             and len(self_crossings_T2([x['D'][0], x['D'][1]])) == 1)
    check(f"each of the 4 resolutions contributes exactly 1 P connector-X (got {nX}/4)",
          nX == 4)
    check("resolved blue still closes up (all components closed)",
          all(is_closed(p) for p in bres))

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
