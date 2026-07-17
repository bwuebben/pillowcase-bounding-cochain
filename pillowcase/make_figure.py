#!/usr/bin/env python3
"""
make_figure.py -- generate paper2's Figure 1 as a PDF via matplotlib, from the
computed curves themselves (RESEARCH_LOG sec 38). Output: ../paper2/fig_q5.pdf,
\\includegraphics'd by main.tex.

The q=5 Lagrangians are folded to the pillowcase fundamental domain
[0,pi] x [0,2pi): blue = R_t(Q_{1/3}+Q_{1/5}) (RAW resolved polyline, 7601 pts,
for smooth rendering), red = R^nat(hat Q_{-1/2}) (both lifts of the earring
figure-eight). Polylines are densified on T^2 and split at genuine fold/wrap
events. Marked: the nine generators (red-blue intersections) and the support
crossings s_A, s_B of the bounding cochain (circled, on blue self-crossings).

NOTE: matplotlib is needed only for THIS figure; all mathematical results in the
repository remain pure-stdlib.
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tangles import TAU, PI, curve, segments, tangle_sum, fiber_circles, south_twists
from resolve import resolve
from polygons import self_intersections_detailed
from earring import f8, P_point
from bigons import simplify, _tdist
from deform import build_geometry, bigon_matrix


def densify(poly, step=0.05):
    out = [poly[0]]
    for a, b in zip(poly, poly[1:]):
        dg = (b[0] - a[0] + PI) % TAU - PI
        dt = (b[1] - a[1] + PI) % TAU - PI
        n = max(1, int(math.ceil(max(abs(dg), abs(dt)) / step)))
        for s in range(1, n + 1):
            out.append(((a[0] + dg * s / n) % TAU, (a[1] + dt * s / n) % TAU))
    return out


def fold_pieces(poly):
    """Fold a T^2 polyline to the fundamental domain, split at fold/wrap events."""
    pts = [P_point(p) for p in densify(poly)]
    pieces, cur = [], [pts[0]]
    for a, b in zip(pts, pts[1:]):
        if abs(b[0] - a[0]) > 0.3 or abs(b[1] - a[1]) > 0.3:
            if len(cur) > 1:
                pieces.append(cur)
            cur = [b]
        else:
            cur.append(b)
    if len(cur) > 1:
        pieces.append(cur)
    return pieces


def main():
    # raw blue (unsimplified) for smooth rendering
    s3 = segments(curve(south_twists(3)))[0]
    s5 = segments(curve(south_twists(5)))[0]
    blue_raw, _ = resolve(tangle_sum(s3, s5), fiber_circles(s3, s5), eps=0.05)
    blue_raw = blue_raw[0]
    red_lifts = f8((2, 1), eps=0.10, phi=0.25)

    # generators and support crossings from the validated pipeline
    red_c, blue_c, _ = build_geometry()
    gens, _d = bigon_matrix(red_c, blue_c)
    gpts = []
    for g in gens:
        pp = P_point(g['pt'])
        if all(_tdist(pp, q) > 1e-3 for q in gpts):
            gpts.append(pp)
    sx = [P_point(s['pt']) for s in self_intersections_detailed(blue_c)]
    sA = min(sx, key=lambda p: _tdist(p, (0.028, 1.272)))
    sB = min(sx, key=lambda p: _tdist(p, (3.057, 4.981)))

    fig, ax = plt.subplots(figsize=(4.6, 8.4))
    # curves
    for piece in fold_pieces(blue_raw):
        xs, ys = zip(*piece)
        ax.plot(xs, ys, color="#1040c0", lw=0.9, solid_joinstyle="round", zorder=2)
    for lift in red_lifts:
        for piece in fold_pieces(lift):
            xs, ys = zip(*piece)
            ax.plot(xs, ys, color="#c01818", lw=0.8, solid_joinstyle="round", zorder=3)
    # pillowcase corners (orbifold points): (0,0),(0,pi),(0,2pi),(pi,0),(pi,pi),(pi,2pi)
    for cg in (0.0, PI):
        for ct in (0.0, PI, TAU):
            ax.plot([cg], [ct], marker="s", ms=4, color="black", zorder=5)
    # generators
    ax.plot([p[0] for p in gpts], [p[1] for p in gpts], "o", ms=4.5,
            color="black", zorder=6)
    # support crossings
    for (p, name, dx, ha) in ((sA, r"$s_A$", 0.10, "left"), (sB, r"$s_B$", -0.10, "right")):
        ax.plot([p[0]], [p[1]], "o", ms=11, mfc="none", mec="black", mew=1.3, zorder=6)
        ax.annotate(name, p, xytext=(p[0] + dx, p[1]), ha=ha, va="center", fontsize=11)
    # frame and labels
    ax.set_xlim(-0.06, PI + 0.06)
    ax.set_ylim(-0.10, TAU + 0.10)
    ax.set_aspect("equal")
    ax.set_xticks([0, PI / 2, PI])
    ax.set_xticklabels([r"$0$", r"$\pi/2$", r"$\pi$"], fontsize=10)
    ax.set_yticks([0, PI, TAU])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"], fontsize=10)
    ax.set_xlabel(r"$\gamma$", fontsize=12)
    ax.set_ylabel(r"$\theta$", fontsize=12, rotation=0, labelpad=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(True)
    ax.tick_params(direction="out", length=3)
    fig.tight_layout(pad=0.4)
    out = "../paper2/fig_q5.pdf"
    fig.savefig(out)
    print(f"wrote {out}: {len(gpts)} generators, "
          f"sA={tuple(round(v,3) for v in sA)}, sB={tuple(round(v,3) for v in sB)}")


if __name__ == "__main__":
    main()
