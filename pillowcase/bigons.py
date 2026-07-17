#!/usr/bin/env python3
"""
bigons.py -- step (d): the pillowcase-aware immersed-bigon counter
(RESEARCH_LOG sec 27).

Combinatorial Floer theory for immersed curves in the pillowcase, de
Silva-Robbin-Salamon style: a Floer bigon (smooth lune) from x to y with
boundary on (red, blue) corresponds to a pair of arcs alpha (on red, x->y) and
beta (on blue, y->x) whose closed loop gamma = alpha . beta has winding function
w with:

  * [gamma] = 0 in H_1(T^2)              (else it bounds nothing),
  * w == 0 at the four iota-fixed points (the pillowcase corners: the disk must
    avoid them, which also anchors the additive constant),
  * w >= 0 on every face (checked on both sides of every edge of gamma),
  * convex corner patterns (m+1, m, m, m) in the four quadrants at x and at y.

Lifting: an immersed bigon in P* lifts to T^2 (simply-connected domain, free
cover away from corners); fixing the lift with alpha on C_+ makes the T^2 count
EQUAL the P count. Blue's corner passages are Lagrangian arc-endpoints: bigon
boundaries may not contain them (this is what kills the unlink's fake (x,y)
lens -- that region is the (x,b,y) mu^2-triangle of the bounding cochain).

Counting convention: an ordered pair (x,y) is accepted when w >= 0 (not all 0);
the reversed loop has w <= 0, so each geometric lune is counted for exactly one
ordering. The differential d[x][y] = #lunes mod 2; rank HF = #gens - 2 rank(d).
Limitation (noted): boundary arcs are single-traversal (no multiply-wrapping
lunes); validated against Smith's unlink (0 bigons) and thm:main (2 bigons).
"""
import math
from tangles import TAU, PI, curve, segments, tangle_sum, fiber_circles, south_twists
from resolve import resolve, _short_edges, is_closed
from earring import f8, arch

CORNERS = [(0.0, 0.0), (0.0, PI), (PI, 0.0), (PI, PI)]


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
def _tdist(p, q):
    dg = abs(p[0] % TAU - q[0] % TAU)
    dt = abs(p[1] % TAU - q[1] % TAU)
    return math.hypot(min(dg, TAU - dg), min(dt, TAU - dt))


def unwrap(poly):
    """Closed T^2 polyline -> continuous plane polyline (cumulative short deltas)."""
    out = [(poly[0][0] % TAU, poly[0][1] % TAU)]
    for p, q in zip(poly, poly[1:]):
        dg = (q[0] - p[0] + PI) % TAU - PI
        dt = (q[1] - p[1] + PI) % TAU - PI
        out.append((out[-1][0] + dg, out[-1][1] + dt))
    return out


def simplify(poly, tol=5e-4):
    """Douglas-Peucker on the unwrapped polyline (keeps torus geometry intact)."""
    u = unwrap(poly)

    def dp(i, j, keep):
        if j <= i + 1:
            return
        ax, ay = u[i]
        bx, by = u[j]
        L = math.hypot(bx - ax, by - ay) or 1e-12
        best, bi = -1.0, None
        for k in range(i + 1, j):
            d = abs((bx - ax) * (u[k][1] - ay) - (by - ay) * (u[k][0] - ax)) / L
            if d > best:
                best, bi = d, k
        if best > tol:
            keep.add(bi)
            dp(i, bi, keep)
            dp(bi, j, keep)

    keep = {0, len(u) - 1}
    import sys
    sys.setrecursionlimit(10000)
    dp(0, len(u) - 1, keep)
    idx = sorted(keep)
    # subdivide long edges: mod-2pi short-edge arithmetic needs |delta| << pi
    out = []
    for a, b in zip(idx, idx[1:]):
        pa, pb = u[a], u[b]
        steps = max(1, int(math.ceil(max(abs(pb[0] - pa[0]),
                                         abs(pb[1] - pa[1])) / 0.8)))
        for s in range(steps):
            lam = s / steps
            out.append(((pa[0] + lam * (pb[0] - pa[0])) % TAU,
                        (pa[1] + lam * (pb[1] - pa[1])) % TAU))
    out.append((u[idx[-1]][0] % TAU, u[idx[-1]][1] % TAU))
    return out


def _seg_cross(a1, b1, a2, b2):
    """Torus-aware transverse segment intersection; returns (t, u, point) or None.
    Segments given in locally-unwrapped short form."""
    mx, my = (a1[0] + b1[0]) / 2, (a1[1] + b1[1]) / 2
    bx, by = (a2[0] + b2[0]) / 2, (a2[1] + b2[1]) / 2
    sx = TAU * round((mx - bx) / TAU)
    sy = TAU * round((my - by) / TAU)
    x1, y1 = a1; x2, y2 = b1
    x3, y3 = a2[0] + sx, a2[1] + sy
    x4, y4 = b2[0] + sx, b2[1] + sy
    d = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3)
    if abs(d) < 1e-14:
        return None
    t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / d
    u = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / d
    if 1e-9 < t < 1 - 1e-9 and 1e-9 < u < 1 - 1e-9:
        return (t, u, (x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
    return None


def edges_of(poly):
    return _short_edges(poly)  # list of (k, a, b) with b locally unwrapped


def intersections_detailed(polyA, polyB):
    """All transverse intersections of two polylines; returns dicts with edge
    indices and parameters on each."""
    EA, EB = edges_of(polyA), edges_of(polyB)
    out = []
    for (k1, a1, b1) in EA:
        for (k2, a2, b2) in EB:
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


def arc_of(poly, kA, tA, kB, tB, forward=True):
    """Polyline arc on the closed polyline from point (kA,tA) to (kB,tB), going
    in the direction of increasing (forward) or decreasing index."""
    n = len(poly) - 1  # closed: last point duplicates first
    E = edges_of(poly)

    def pt_on(k, t):
        _, a, b = E[k]
        return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))

    pts = [pt_on(kA, tA)]
    if forward:
        k = kA
        while True:
            k2 = (k + 1) % n
            pts.append(poly[k2])
            if k2 == kB and True:
                pass
            k = k2
            if k == kB:
                break
            if k == kA:  # wrapped fully without finding kB (kB==kA handled below)
                break
        pts.append(pt_on(kB, tB))
        if kA == kB and tB > tA:
            pts = [pt_on(kA, tA), pt_on(kB, tB)]
    else:
        k = kA
        while True:
            pts.append(poly[k])
            k = (k - 1) % n
            if (k + 1) % n == (kB + 1) % n and k == kB:
                break
            if k == kA:
                break
        pts.append(pt_on(kB, tB))
        if kA == kB and tB < tA:
            pts = [pt_on(kA, tA), pt_on(kB, tB)]
    return pts


def contains_vertex(poly, kA, kB, forward, special_ks):
    """Does the forward/backward arc from edge kA to edge kB pass through any
    vertex index in special_ks (vertex j = start of edge j)?"""
    n = len(poly) - 1
    ks = set(special_ks)
    if forward:
        k = (kA + 1) % n
        while k != (kB + 1) % n:
            if k in ks:
                return True
            k = (k + 1) % n
            if k == (kA + 1) % n:
                break
    else:
        k = kA
        while k != kB:
            if k in ks:
                return True
            k = (k - 1) % n
            if k == kA:
                break
    return False


# ---------------------------------------------------------------------------
# winding machinery
# ---------------------------------------------------------------------------
def loop_edges(loop_pts):
    out = []
    for p, q in zip(loop_pts, loop_pts[1:]):
        dg = (q[0] - p[0] + PI) % TAU - PI
        dt = (q[1] - p[1] + PI) % TAU - PI
        if abs(dg) + abs(dt) > 1e-12:
            g0, t0 = p[0] % TAU, p[1] % TAU
            out.append(((g0, t0), (g0 + dg, t0 + dt)))
    return out


def loop_class(loop_pts):
    """Homology class of the loop in H_1(T^2) (must be (0,0) to bound)."""
    G = T = 0.0
    for p, q in zip(loop_pts, loop_pts[1:]):
        G += (q[0] - p[0] + PI) % TAU - PI
        T += (q[1] - p[1] + PI) % TAU - PI
    return (round(G / TAU), round(T / TAU), abs(G / TAU - round(G / TAU)) +
            abs(T / TAU - round(T / TAU)))


def winding_at(z, E, zref, jitters=((0.0, 0.0), (1.3e-4, 0.7e-4), (-0.9e-4, 1.1e-4),
                                    (0.5e-4, -1.4e-4))):
    """w(z) - w(zref) by signed crossings along the short segment z -> zref
    (valid because [gamma]=0). Retries with jittered z on near-degenerate hits."""
    for (jx, jy) in jitters:
        zz = (z[0] + jx, z[1] + jy)
        dx = (zref[0] - zz[0] + PI) % TAU - PI
        dy = (zref[1] - zz[1] + PI) % TAU - PI
        a1, b1 = (zz[0] % TAU, zz[1] % TAU), (zz[0] % TAU + dx, zz[1] % TAU + dy)
        tot, bad = 0, False
        for (a2, b2) in E:
            mx, my = (a1[0] + b1[0]) / 2, (a1[1] + b1[1]) / 2
            bx, by = (a2[0] + b2[0]) / 2, (a2[1] + b2[1]) / 2
            sx = TAU * round((mx - bx) / TAU)
            sy = TAU * round((my - by) / TAU)
            x1, y1 = a1; x2, y2 = b1
            x3, y3 = a2[0] + sx, a2[1] + sy
            x4, y4 = b2[0] + sx, b2[1] + sy
            d = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3)
            if abs(d) < 1e-14:
                continue
            t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / d
            u = ((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / d
            if -1e-7 < t < 1e-7 or 1 - 1e-7 < t < 1 + 1e-7 or \
               -1e-7 < u < 1e-7 or 1 - 1e-7 < u < 1 + 1e-7:
                bad = True
                break
            if 0 < t < 1 and 0 < u < 1:
                tot += 1 if (x4 - x3) * (y2 - y1) - (y4 - y3) * (x2 - x1) > 0 else -1
        if not bad:
            return tot
    return tot  # last attempt


def _pt_seg_dist(p, a, b):
    vx, vy = b[0] - a[0], b[1] - a[1]
    wx, wy = p[0] - a[0], p[1] - a[1]
    L2 = vx * vx + vy * vy or 1e-24
    t = max(0.0, min(1.0, (wx * vx + wy * vy) / L2))
    return math.hypot(wx - t * vx, wy - t * vy)


def _corner_ok(E, c0, p, u_dir, v_dir, r_clear):
    """Convex-corner test at loop corner p with outgoing rays u,v (both FROM p):
    the two sectors' windings must be {m, m+1}, m >= 0, with the (m+1)-side (the
    disk) subtending angle < pi. Adaptive sampling keeps clearance from the rays."""
    au = math.atan2(u_dir[1], u_dir[0])
    av = math.atan2(v_dir[1], v_dir[0])
    rel = (av - au) % TAU                      # sector 1: CCW from u to v
    vals, opens = [], []
    for (start, opening) in ((au, rel), (av, TAU - rel)):
        if opening < 1e-4:
            return False                        # degenerate tangency
        sig = opening / 2
        delta = max(1.2e-3, 3.2e-4 / max(math.sin(min(sig, PI / 2)), 1e-3))
        if r_clear > 0:
            delta = min(delta, 0.45 * r_clear) or delta
        delta = max(delta, 3e-4)
        d = (math.cos(start + sig), math.sin(start + sig))
        vals.append(winding_at((p[0] + delta * d[0], p[1] + delta * d[1]), E, c0))
        opens.append(opening)
    lo, hi = min(vals), max(vals)
    if not (hi == lo + 1 and lo >= 0):
        return False
    disk_opening = opens[vals.index(hi)]
    return disk_opening < PI


def is_lune(alpha, beta_rev, x, y, dirs_x=None, dirs_y=None):
    """Check whether the loop alpha (x->y on red) + beta_rev (y->x on blue)
    bounds an immersed bigon (smooth lune with convex corners)."""
    loop = list(alpha) + list(beta_rev)[1:]
    g, t, frac = loop_class(loop)
    if g != 0 or t != 0 or frac > 1e-3:
        return False
    E = loop_edges(loop)
    if len(E) < 3:
        return False
    # anchor: w(corner_0) = 0; all pillowcase corners must agree
    c0 = CORNERS[0]
    for c in CORNERS[1:]:
        if winding_at(c, E, c0) != 0:
            return False
    # per-edge clearance (nearest non-adjacent edge) for adaptive side sampling
    nE = len(E)
    mids = [((a[0] + b[0]) / 2, (a[1] + b[1]) / 2) for (a, b) in E]
    wmax = 0
    for i, (a, b) in enumerate(E):
        clear = 1.0
        for j, (a2, b2) in enumerate(E):
            dj = min((i - j) % nE, (j - i) % nE)
            if dj <= 1:
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
    # convex corners at x (loop start/end) and y (the alpha->beta junction)
    kx = 0
    ky = len(list(alpha)) - 1
    lp = loop

    def ray(i, sign):
        p = lp[i]
        q = lp[(i + sign) % (len(lp) - 1)] if sign > 0 else lp[i - 1]
        dg = (q[0] - p[0] + PI) % TAU - PI
        dt = (q[1] - p[1] + PI) % TAU - PI
        return (dg, dt)

    def clear_at(p, iskip):
        c = 1.0
        for j, (a2, b2) in enumerate(E):
            if min((iskip - j) % nE, (j - iskip) % nE) <= 1:
                continue
            c = min(c, _pt_seg_dist(p, a2, b2))
        return c

    ux = ray(0, +1)                             # out along alpha
    vx = ray(len(lp) - 1, -1)                   # out along (reversed) end of loop
    if not _corner_ok(E, c0, (lp[0][0] % TAU, lp[0][1] % TAU), ux, vx,
                      clear_at((lp[0][0] % TAU, lp[0][1] % TAU), 0)):
        return False
    uy = ray(ky, -1)                            # back along alpha
    vy = ray(ky, +1)                            # out along beta_rev
    if not _corner_ok(E, c0, (lp[ky][0] % TAU, lp[ky][1] % TAU), uy, vy,
                      clear_at((lp[ky][0] % TAU, lp[ky][1] % TAU), ky)):
        return False
    return True


# ---------------------------------------------------------------------------
# the Floer differential between two curves
# ---------------------------------------------------------------------------
def floer_rank(redC, blue, verbose=False):
    """redC: ONE lift (C_+) of the red immersed P-curve (closed polyline);
    blue: the blue T^2 curve (closed polyline, iota-invariant). Counts lunes
    between all generator pairs; returns (n_gens, n_lunes, rank_d, rank_HF)."""
    hits = intersections_detailed(redC, blue)
    n = len(hits)
    # blue corner vertices (Lagrangian arc endpoints)
    ncv = [k for k, p in enumerate(blue[:-1])
           if any(_tdist(p, c) < 1e-6 for c in CORNERS)]
    lunes = []
    d = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            x, y = hits[i], hits[j]
            for fa in (True, False):
                alpha = arc_of(redC, x['kA'], x['tA'], y['kA'], y['tA'], fa)
                for fb in (True, False):
                    if contains_vertex(blue, y['kB'], x['kB'], fb, ncv):
                        continue
                    beta_rev = arc_of(blue, y['kB'], y['tB'], x['kB'], x['tB'], fb)
                    if is_lune(alpha, beta_rev, x['pt'], y['pt'],
                               (x['dirA'], x['dirB']), (y['dirA'], y['dirB'])):
                        d[i][j] ^= 1
                        lunes.append((i, j, fa, fb))
                        if verbose:
                            print(f"      lune: gen{i} {x['pt']} -> gen{j} {y['pt']}"
                                  f"  (red fwd={fa}, blue fwd={fb})")
    # rank over F2
    M = [row[:] for row in d]
    rank = 0
    for col in range(n):
        piv = next((r for r in range(rank, n) if M[r][col]), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        for r in range(n):
            if r != rank and M[r][col]:
                M[r] = [a ^ b for a, b in zip(M[r], M[rank])]
        rank += 1
    return n, len(lunes), rank, n - 2 * rank


# ---------------------------------------------------------------------------
# battery (sec 27)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    results = []

    def check(name, ok):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    print("== the pillowcase Floer bigon counter (RESEARCH_LOG sec 27) ==\n")

    print("(16) synthetic: two circles crossing twice, corners all outside")
    NA = 120
    cA = [(2.0 + 0.6 * math.cos(TAU * k / NA), 2.2 + 0.6 * math.sin(TAU * k / NA))
          for k in range(NA + 1)]
    cB = [(2.8 + 0.59 * math.cos(TAU * k / NA + 0.13),
           2.217 + 0.59 * math.sin(TAU * k / NA + 0.13)) for k in range(NA + 1)]
    n, nl, r, hf = floer_rank(cA, cB)
    check(f"2 gens, 3 lunes (lens + 2 crescents), d(x)=y, d(y)=2x=0, HF = 0 "
          f"(got {n} gens, {nl} lunes, rank {r}, HF {hf})",
          n == 2 and nl == 3 and r == 1 and hf == 0)

    print("\n(17) the unlink (Smith Fig UnlinkA): NO bigons -- the (x,b,y) region"
          " is a mu^2-triangle, not a bigon")
    red_u = f8((1, 0), eps=0.30, phi=0.20, N=601)
    blue_u = arch(amp=0.5)[0]
    ru = simplify(red_u[0], 2e-4)
    bu = simplify(blue_u, 2e-4)
    n, nl, r, hf = floer_rank(ru, bu)
    check(f"2 generators, 0 bigons, rank HF = 2 (got {n} gens, {nl} lunes, "
          f"HF {hf})", n == 2 and nl == 0 and hf == 2)

    print("\n(18) THE GATE: P(-2,3,5) = red x blue -> 2 bigons, rank HF = 5")
    s3 = segments(curve(south_twists(3)))[0]
    s5 = segments(curve(south_twists(5)))[0]
    blue_polys, _ = resolve(tangle_sum(s3, s5), fiber_circles(s3, s5), eps=0.05)
    assert len(blue_polys) == 1
    blue = simplify(blue_polys[0], 3e-4)
    red = simplify(f8((2, 1), eps=0.10, phi=0.25)[0], 2e-4)
    n, nl, r, hf = floer_rank(red, blue, verbose=True)
    print(f"        generators {n}, lunes {nl}, rank(d) {r}, rank HF {hf}")
    check(f"9 generators on the C+ lift (got {n})", n == 9)
    check(f"2 bigons (Smith thm:main; got {nl})", nl == 2)
    check(f"rank HF = 5 (got {hf})", hf == 5)

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
