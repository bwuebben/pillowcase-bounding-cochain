#!/usr/bin/env python3
"""
earring.py -- step (c): the earring (figure-eight) curve of a rational tangle
(RESEARCH_LOG sec 26).

Decoded from Smith's Unlink figures + HK's model formula, now understood:
R^natural of a rational tangle whose pillowcase curve is the line ell(t) = t*(q,p)
through two corners is the DOUBLED arc with ONE swap:

    C_+ = { ell(t) + eps*cos(t+phi)*nhat : t in [0,2pi) },   C_- = iota(C_+),

where nhat is the unit normal. The two copies' ends connect through the pillowcase
fold edges AROUND the two corners (the fold identification (0,eps)~(0,2pi-eps) does
the corner-looping for free -- this is why HK's formula [sigma, -2s cos sigma]
needed no extra symplectomorphism: the sec-20 "Psi wall" was a misreading of the
fold). The single P self-crossing (the F8 pinch, support of the would-be b, which
is 0 for rational tangles by Smith lem:trivialbc) is FORCED to the arc's parameter
midpoint t = pi/2; phi shifts it slightly off the base line (keep sin(phi) != 0).

Pairing frame (derived in sec 26 and confirmed by Smith's hat: phi^* intertwines
with diag(1,-1), his line 763): for K = num(T1+T2) split along T1's Conway sphere,
the outside tangle's curve enters the pairing MIRRORED, (gamma,theta) -> (gamma,
-theta). Equivalently (Smith's choice, adopted here): red = F8 of the mirrored
base line of Q_{-1/2} = the slope +1/2 line, blue = R_t(Q_{1/3}+Q_{1/5}) as built.
Det check: |2*8 - 1*15| = 1 = det P(-2,3,5), one corner point (the abelian rep).
"""
import math
from tangles import TAU, PI, curve, segments, tangle_sum, fiber_circles, south_twists
from resolve import (assemble, resolve, self_crossings_T2, fold_crossings,
                     _seg_int_torus, _short_edges, is_closed)


def f8(qp, eps=0.10, phi=0.25, N=1201, jit=0.171):
    """The earring figure-eight of the corner-to-corner tangle line t -> t*(q,p):
    T^2 representative = the pair {C_+, iota(C_+)}. The sample jitter keeps the
    forced pinch (t = pi/2) off the sample grid (vertex-coincident crossings are
    invisible to the interior-only intersection test)."""
    q, p = qp
    L = math.hypot(q, p)
    ng, nt = -p / L, q / L
    Cp = []
    for k in range(N + 1):
        t = TAU * k / N + jit
        w = eps * math.cos(t + phi)
        Cp.append(((q * t + w * ng) % TAU, (p * t + w * nt) % TAU))
    Cm = [((-g) % TAU, (-t_) % TAU) for (g, t_) in Cp]
    return [Cp, Cm]


def fold_pairs(pts, tol=1e-3):
    """Exact iota-pairing of T^2 crossing points: each point is matched with the
    nearest iota-image among the rest. Returns (n_pairs, worst_match_distance)."""
    def d(p, q):
        dg = abs(p[0] % TAU - q[0] % TAU)
        dt = abs(p[1] % TAU - q[1] % TAU)
        return math.hypot(min(dg, TAU - dg), min(dt, TAU - dt))
    left = list(pts)
    pairs, worst = 0, 0.0
    while left:
        p = left.pop()
        ip = ((-p[0]) % TAU, (-p[1]) % TAU)
        j, best = min(enumerate(d(ip, q) for q in left), key=lambda x: x[1],
                      default=(None, None))
        if j is None:
            return (pairs, float('inf'))
        worst = max(worst, best)
        left.pop(j)
        pairs += 1
    return (pairs, worst)


def P_point(p):
    """Canonical pillowcase representative of a T^2 point."""
    g, t = p[0] % TAU, p[1] % TAU
    if g > PI + 1e-12 or (abs(g) < 1e-12 and t > PI):
        g, t = (TAU - g) % TAU, (TAU - t) % TAU
    return (g, t)


def arch(amp=0.35, N=400):
    """Sheared R(Q_0): the iota-invariant arch theta = amp*sin(gamma) (Smith
    Fig UnlinkA blue)."""
    return [[(TAU * k / N, (amp * math.sin(TAU * k / N)) % TAU) for k in range(N + 1)]]


def curve_intersections(polysA, polysB):
    """Transverse intersections between two curve families on T^2."""
    EA = [(pi_, k, a, b) for pi_, poly in enumerate(polysA)
          for (k, a, b) in _short_edges(poly)]
    EB = [(pi_, k, a, b) for pi_, poly in enumerate(polysB)
          for (k, a, b) in _short_edges(poly)]
    pts = []
    for (_, _, a1, b1) in EA:
        for (_, _, a2, b2) in EB:
            r = _seg_int_torus(a1, b1, a2, b2)
            if r:
                pts.append((r[0] % TAU, r[1] % TAU))
    return pts


if __name__ == "__main__":
    results = []

    def check(name, ok):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    print("== the earring curve (RESEARCH_LOG sec 26) ==\n")

    print("(14) calibration: the unlink U_2 = Q_0 u Q_0^natural  [Smith Fig UnlinkA]")
    red_u = f8((1, 0), eps=0.30, phi=0.20)
    xr = self_crossings_T2(red_u)
    nP, worst = fold_pairs(xr)
    check(f"R^nat(Q_0) has exactly 1 P self-crossing (got {nP}, pair-dist {worst:.1e})",
          nP == 1 and worst < 1e-3)
    if xr:
        pp = P_point(xr[0])
        tp = min(pp[1], TAU - pp[1])
        check(f"pinch at the arc midpoint (pi/2, 0)-ish (got ({pp[0]:.3f},{tp:.3f}))",
              abs(pp[0] - PI / 2) < 0.05 and tp < 0.12)
    blue_u = arch(amp=0.5)
    xi = curve_intersections(red_u, blue_u)
    nPi, worst = fold_pairs(xi)
    check(f"|R^nat(Q_0) ^ sheared R(Q_0)| = 2 P-points x,y (got {nPi}, T^2 {len(xi)})",
          nPi == 2 and len(xi) == 4)

    print("\n(15) THE GATE: red x blue for P(-2,3,5)  [Smith thm:main: 9 points]")
    # blue as built in resolve.py
    s3 = segments(curve(south_twists(3)))[0]
    s5 = segments(curve(south_twists(5)))[0]
    bsum = tangle_sum(s3, s5)
    bcirc = fiber_circles(s3, s5)
    blue, _ = resolve(bsum, bcirc, eps=0.05)
    # red: F8 of the mirrored Q_{-1/2} line = direction (2,+1), Smith's hat frame
    red = f8((2, 1), eps=0.10, phi=0.25)
    xrb = curve_intersections(red, blue)
    nP, worst = fold_pairs(xrb)
    print(f"        T^2 intersections: {len(xrb)}; P intersections: {nP} "
          f"(iota pair-dist {worst:.1e})")
    check("iota-pairing exact (even count, partners within 1e-3)",
          len(xrb) % 2 == 0 and worst < 1e-3)
    check(f"|red ^ blue| = 9 (Smith)  (got {nP})", nP == 9)
    seen = []
    for p in xrb:
        pp = P_point(p)
        if all(math.hypot(pp[0] - q[0], min(abs(pp[1] - q[1]),
               TAU - abs(pp[1] - q[1]))) > 2e-3 for q in seen):
            seen.append(pp)
            print(f"        generator at P(gamma,theta) = ({pp[0]:.4f}, {pp[1]:.4f})")

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
