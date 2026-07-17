#!/usr/bin/env python3
"""
deform.py -- step (e), part 2: the deformed differential of P(-2,3,5)
(RESEARCH_LOG sec 28).

For each blue self-crossing s and each ordered generator pair (x,y), count the
immersed triangles that contribute to the deformed differential partial_b:

    partial_b(x) = sum_y [ #bigons(x,y) + sum_i c_i #tri(x,s_i,y)
                           + (higher polygons) ] y   (mod 2),   b = sum_i c_i s_i.

A triangle for partial_b(x)->y has corners x, y (generators = red^blue) and s (a
blue self-crossing); its boundary is  red arc(x->y) + blue arc(y->s) + blue arc
(s->x), switching blue branch at s. Counting on the C_+ lift with the red arc on
C_+ gives the pillowcase count directly (same convention validated for the bigons
in sec 27). rank(HF) = #gens - 2 rank(partial_b); target 7 needs rank 1, i.e. the
b-triangles must cancel exactly one of the two rank-5 bigons.

This module builds, per self-crossing s, the F_2 matrix Tri[s][i][j] = #tri(g_i ->
g_j through s) mod 2, and reports which crossings flip a bigon entry.
"""
import math
from tangles import TAU, PI, curve, segments, tangle_sum, fiber_circles, south_twists
from resolve import resolve
from earring import f8, P_point
from bigons import (simplify, intersections_detailed, arc_of, is_lune,
                    contains_vertex, floer_rank, _tdist, CORNERS)
from polygons import self_intersections_detailed, arc_between, _assemble_loop, bounds_disk


def build_geometry():
    return build_geometry_p()


def build_geometry_p(blue_eps=0.05, red_eps=0.10, red_phi=0.25):
    """Parametric geometry -- for perturbation-stability checks (RESEARCH_LOG sec 30)."""
    return build_pretzel(2, blue_eps, red_eps, red_phi)


def build_pretzel(k, blue_eps=0.05, red_eps=0.16, red_phi=0.40):
    """The pillowcase Lagrangians for P(-2,3,2k+1) = num(Q_{-1/2}+Q_{1/3}+Q_{1/(2k+1)})
    (RESEARCH_LOG sec 33). blue = R_t(Q_{1/3}+Q_{1/(2k+1)}) (the perturbed Conway sum,
    resolved); red = R^natural(hat Q_{-1/2}) (the earring, FIXED across the family).
    Requires gcd(2k+1,3)=1 (else corner circles appear and resolve fails). NOTE the
    default red perturbation (0.16,0.40) is GENERIC for the family; (0.10,0.25) is
    degenerate for k=3 (collapses 13 generators to 9)."""
    q = 2 * k + 1
    assert (q % 3) != 0, f"P(-2,3,{q}) needs gcd({q},3)=1 (corner-circle-free)"
    s3 = segments(curve(south_twists(3)))[0]
    sq = segments(curve(south_twists(q)))[0]
    circ = fiber_circles(s3, sq)
    blue_polys, xinfo = resolve(tangle_sum(s3, sq), circ, eps=blue_eps)
    assert len(blue_polys) == 1, f"blue has {len(blue_polys)} components"
    blue = simplify(blue_polys[0], 3e-4)
    red = simplify(f8((2, 1), eps=red_eps, phi=red_phi)[0], 2e-4)
    return red, blue, xinfo


def bigon_matrix(red, blue):
    """Return (gens, d) where gens = list of generator hit-dicts and d[i][j] =
    #bigons(g_i -> g_j) mod 2 (the undeformed differential)."""
    hits = intersections_detailed(red, blue)
    n = len(hits)
    ncv = [k for k, p in enumerate(blue[:-1])
           if any(_tdist(p, c) < 1e-6 for c in CORNERS)]
    d = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            x, y = hits[i], hits[j]
            for fa in (True, False):
                alpha = arc_of(red, x['kA'], x['tA'], y['kA'], y['tA'], fa)
                for fb in (True, False):
                    if contains_vertex(blue, y['kB'], x['kB'], fb, ncv):
                        continue
                    beta = arc_of(blue, y['kB'], y['tB'], x['kB'], x['tB'], fb)
                    if is_lune(alpha, beta, x['pt'], y['pt']):
                        d[i][j] ^= 1
    return hits, d


def _arc_short(poly, k0, k1, forward, maxspan):
    """True if the directed arc on a closed poly from edge k0 to edge k1 spans at
    most maxspan edges (prunes globe-wrapping arcs -- immersed triangles are local)."""
    n = len(poly) - 1
    span = (k1 - k0) % n if forward else (k0 - k1) % n
    return span <= maxspan


def triangle_contributions(red, blue, gens, maxspan_blue=140, maxspan_red=60,
                           only_s=None, verbose=False):
    """For every blue self-crossing s and ordered generator pair (i,j), count
    mod-2 the immersed triangles g_i -> g_j through s. Returns (scross, Tri) where
    scross = list of self-crossing dicts (P-orbit representatives) and Tri is a
    dict s_index -> n x n F_2 matrix."""
    scross_all = self_intersections_detailed(blue)
    # fold to P-orbit representatives (each P self-crossing has 2 T^2 preimages)
    reps = []
    for s in scross_all:
        pp = P_point(s['pt'])
        if only_s is not None and not any(_tdist(pp, q) < 3e-3 for q in only_s):
            continue
        reps.append(s)
    n = len(gens)
    Tri = {}
    for si, s in enumerate(reps):
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                x, y = gens[i], gens[j]
                cnt = 0
                for (kIn, tIn, kOut, tOut) in (
                        (s['kA'], s['tA'], s['kB'], s['tB']),
                        (s['kB'], s['tB'], s['kA'], s['tA'])):
                    for fr in (True, False):
                        if not _arc_short(red, x['kA'], y['kA'], fr, maxspan_red):
                            continue
                        red_arc = arc_between(red, x['kA'], x['tA'],
                                              y['kA'], y['tA'], fr)
                        for f1 in (True, False):
                            if not _arc_short(blue, y['kB'], kIn, f1, maxspan_blue):
                                continue
                            blue1 = arc_between(blue, y['kB'], y['tB'], kIn, tIn, f1)
                            for f2 in (True, False):
                                if not _arc_short(blue, kOut, x['kB'], f2, maxspan_blue):
                                    continue
                                blue2 = arc_between(blue, kOut, tOut,
                                                    x['kB'], x['tB'], f2)
                                loop, corners = _assemble_loop([red_arc, blue1, blue2])
                                if bounds_disk(loop, corners):
                                    cnt += 1
                                    if verbose:
                                        print(f"      tri g{i}->g{j} thru s{si} "
                                              f"P{P_point(s['pt'])}")
                M[i][j] = cnt % 2
        Tri[si] = (s, M)
    return reps, Tri


def rank_f2(M):
    n = len(M)
    A = [row[:] for row in M]
    rank = 0
    for col in range(len(A[0]) if A else 0):
        piv = next((r for r in range(rank, n) if A[r][col]), None)
        if piv is None:
            continue
        A[rank], A[piv] = A[piv], A[rank]
        for r in range(n):
            if r != rank and A[r][col]:
                A[r] = [a ^ b for a, b in zip(A[r], A[rank])]
        rank += 1
    return rank


def triangle_contributions_P(red, blue, gens, maxspan_blue=180, maxspan_red=60):
    """Full sweep: for every blue self-crossing (all T^2 preimages) count mod-2
    the triangles g_i -> g_j through it, then GROUP the two T^2 preimages of each
    pillowcase self-crossing into one P-orbit and UNION their triangle entries
    (the P-count with red arc on C_+ picks whichever preimage closes; different
    P-triangles through one P-crossing may use different preimages).
    Returns (Pcross, TriP): Pcross = list of (Prep_point, [preimage dicts]);
    TriP = list of n x n F_2 matrices, one per P-crossing."""
    scross = self_intersections_detailed(blue)
    n = len(gens)
    # per T^2 preimage, its triangle matrix
    mats = []
    for s in scross:
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                x, y = gens[i], gens[j]
                cnt = 0
                for (kIn, tIn, kOut, tOut) in (
                        (s['kA'], s['tA'], s['kB'], s['tB']),
                        (s['kB'], s['tB'], s['kA'], s['tA'])):
                    for fr in (True, False):
                        if not _arc_short(red, x['kA'], y['kA'], fr, maxspan_red):
                            continue
                        red_arc = arc_between(red, x['kA'], x['tA'],
                                              y['kA'], y['tA'], fr)
                        for f1 in (True, False):
                            if not _arc_short(blue, y['kB'], kIn, f1, maxspan_blue):
                                continue
                            blue1 = arc_between(blue, y['kB'], y['tB'], kIn, tIn, f1)
                            for f2 in (True, False):
                                if not _arc_short(blue, kOut, x['kB'], f2, maxspan_blue):
                                    continue
                                blue2 = arc_between(blue, kOut, tOut,
                                                    x['kB'], x['tB'], f2)
                                loop, corners = _assemble_loop([red_arc, blue1, blue2])
                                if bounds_disk(loop, corners):
                                    cnt += 1
                M[i][j] = cnt % 2
        mats.append((s, M))
    # group by P-orbit (iota-partners)
    Pcross, TriP = [], []
    used = [False] * len(mats)
    for a in range(len(mats)):
        if used[a]:
            continue
        sa, Ma = mats[a]
        pa = P_point(sa['pt'])
        group = [Ma]
        preims = [sa]
        used[a] = True
        for b in range(a + 1, len(mats)):
            if used[b]:
                continue
            sb, Mb = mats[b]
            if _tdist(P_point(sb['pt']), pa) < 3e-3:
                group.append(Mb)
                preims.append(sb)
                used[b] = True
        U = [[0] * n for _ in range(n)]
        for M in group:
            for i in range(n):
                for j in range(n):
                    U[i][j] ^= M[i][j]
        Pcross.append((pa, preims))
        TriP.append(U)
    return Pcross, TriP


if __name__ == "__main__":
    import json, sys
    red, blue, xinfo = build_geometry()
    gens, d = bigon_matrix(red, blue)
    n = len(gens)
    print(f"generators: {n}; undeformed rank(d) = {rank_f2(d)}  "
          f"HF = {n - 2 * rank_f2(d)}")
    bigon_pairs = [(i, j) for i in range(n) for j in range(n) if d[i][j]]
    print(f"bigon entries d[i][j]=1: {bigon_pairs}")

    print("\n=== full triangle sweep over all 40 P self-crossings ===")
    Pcross, TriP = triangle_contributions_P(red, blue, gens)
    print(f"P self-crossings: {len(Pcross)}")
    active = []
    for k, (pp, preims) in enumerate(Pcross):
        ent = [(i, j) for i in range(n) for j in range(n) if TriP[k][i][j]]
        if ent:
            active.append(k)
            print(f"  S{k} @P{tuple(round(v,4) for v in pp)} "
                  f"({len(preims)} preim): triangles {ent}")
    print(f"\n{len(active)} of {len(Pcross)} P-crossings carry triangles")

    # save for the solver
    out = dict(n=n, d=d, bigons=bigon_pairs,
               gens=[list(P_point(g['pt'])) for g in gens],
               Pcross=[list(pp) for pp, _ in Pcross],
               TriP=TriP)
    with open("deform_data.json", "w") as f:
        json.dump(out, f)
    print("wrote deform_data.json")
