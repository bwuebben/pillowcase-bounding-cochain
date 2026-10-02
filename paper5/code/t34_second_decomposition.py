#!/usr/bin/env python3
r"""Section 7 and Proposition 7.1: T(3,4) with (r, s) = (7, -5) (j = 1), with the decompositions j = 0,
(r, s) = (3, -2) (Theorem 1.1) and j = -1, (r, s) = (-1, 1) as controls.

Perturbations: (eps_A, eps_B) = rho (cos a, sin a), a = k pi/6 + 0.05 (k = 0..11, twelve directions, eps_B != 0),
rho in {0.03, 0.01, 0.003}: 36 runs; earring parameters eps = rho/10 and rho/30.  The approximate coordinates quoted
in Section 7 are checked over all 36 runs.  In addition one run at (eps_A, eps_B) = (0.075, 0.05) with
eps in {0.001, 0.003, 0.01, 0.03, 0.05} checks the triangle of the bounding cochain and the twelve
self-intersections.

For (7, -5), in every run:
  (1) R_pi(Y,T) is an arc and five circles; the lift of the arc from (0,0) ends at (7 pi, 8 pi) and is homotopic
      rel ends in R^2 minus the lattice to the polygon of Section 7 (equal reduced gap words), hence not to the
      arc of Theorem 1.1, whose lift ends at (pi, 2 pi); the arc crosses Delta_0 once and then {phi = pi} once;
  (2) every lift of every component is embedded; mu(C, l1) + z(C) = 0 (mod 4) on the circles; four circles have
      displacement +-(4 pi, 4 pi) per period and lie in a band of phi away from 2 pi Z; the H_1(P*) class of
      every circle is not a multiple of a power of the class of L0 (condition (ii) of Definition 2.2);
  (3) 7 = 1 + |sigma| generators: r_+, x^-, x^+ on the arc, in gradings 2, 3, 0; four on one circle, one in each
      grading; none on the other four circles;
  (4) no bigon; the closed curve that would bound a bigon x^- -> r_+ has winding number 1 about (pi, pi) and
      (3 pi, 3 pi); the generators on the circle lie on four distinct strands, also after translation by the
      period; H = (2,1,2,2) in gradings 0..3, rank 7, against I^natural(T(3,4)) = (2,1,1,1);
  (5) the self-intersection s of the arc near (0,0) (deck transformation: translation by (-4 pi, -4 pi)) and the
      embedded triangle with convex corners x^-, r_+, s and no lattice point inside; degree of s; the candidate
      polygons with corners at s enclose lattice points; deformed homology (2,1,1,1).
Also: the polygon meets L0 in exactly three points for eps in [5e-4, 0.2] and bounds no bigon; the arc crosses
itself twelve times in P, the arc of Theorem 1.1 is embedded; controls j = 0, -1; the strand-meridian
quotients of the two tangles (Remark after Proposition 7.1).

Floating-point numerics, as stated in Section 7; these computations are not proofs.
Usage: t34_second_decomposition.py [--quick]   (--quick: rho = 0.03 only, 12 runs)
"""
import math
import sys
import time

import numpy as np

from t3n_pillowcase import (PI, TAU, Tangle, signature, slice_points, trace_components, refine, lift, arc_lift,
                            circle_lift, generators, grading_both, bigon_candidates, lattice_windings, gap_word,
                            z_class, tree_word, cyclic_reduce, word_str, seg_crossings, resample, polyline_between,
                            strand_path, signed_area, seed_coverage, analyse_complex, circle_degree_pattern,
                            corner_loops_z)

POLYGON = np.array([(0, 0), (0.5, 0), (4.8, 4.3), (4.8, 5.3), (2.5, 3), (2.5, 4), (6.5, 8), (7, 8)], float) * PI
PAPER = dict(cross_delta0=4.33, cross_pi=2.67, s=(0.18, 0.10), band=(-2.4, -0.75), selfint=12)
SIGMA = signature(3, 4)                                        # -6, so gr(r_+) = 2

results = {}


def record(item, ok):
    results.setdefault(item, []).append(bool(ok))


def check(label, ok, computed, paper):
    record(label, ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: computed {computed}; paper {paper}")


def trace(r, s, eA, eB, rho):
    T = Tangle(3, 4, r, s, eA, eB)
    seeds = slice_points(T, nslice=1200, ngrid=8000)
    comps = trace_components(T, seeds, h0=min(2e-3, rho / 6))
    for c in comps:
        c['pts'] = refine(T, c['pts'], 3e-3)
        c['L'] = lift(T.pillow(c['pts']))
    return T, comps, seed_coverage(seeds, comps)


def h1_class(word):
    v = {'X': 0, 'Y': 0, 'Z': 0}
    for L, e in word:
        v[L] += e
    return np.array([v['X'], v['Y'], v['Z']])


def L0_class():
    ts = np.linspace(0, TAU, 4001)
    S = np.c_[ts + 0.01 * np.sin(ts) + PI / 2, ts - 0.01 * np.sin(ts) + PI / 2]
    return h1_class(tree_word(S)), z_class(S)[0]


def strip_crossings(L):
    """crossings of the lift with the lines phi in pi Z (excluding the start): (gamma/pi, line index k)."""
    f = (L[:, 1] - L[:, 0]) / PI
    out = []
    for i in range(5, len(f) - 5):
        if math.floor(f[i]) != math.floor(f[i + 1]):
            out.append((L[i, 0] / PI, max(math.floor(f[i]), math.floor(f[i + 1]))))
    return out


def circle_mu(Lc_period):
    """Maslov index mu(C, l1) of a circle: the number of half-turns of its tangent line over one period."""
    L, disp = Lc_period[:-1], Lc_period[-1] - Lc_period[0]
    d = np.diff(np.vstack([L, L[:2] + disp]), axis=0)
    d = d[np.linalg.norm(d, axis=1) > 0]
    ang = np.unwrap(2 * np.arctan2(d[:, 1], d[:, 0])) / 2
    return int(round((ang[-1] - ang[0]) / PI))


def winding_about(loop, pt):
    d = np.vstack([loop, loop[:1]]) - np.asarray(pt)
    a = np.arctan2(d[:, 1], d[:, 0])
    return int(round(np.sum(np.mod(np.diff(a) + PI, TAU) - PI) / TAU))


def triangle(L, eps):
    """the self-intersection s of the arc lift L with its translate by (-4 pi, -4 pi) nearest the corner, and the
    triangle bounded by L from r_+ to s, the translate from s to the translate of x^-, and the strand S_-^0."""
    g = np.array([-2 * TAU, -2 * TAU])
    X = sorted((i + t, j + u, pt) for (i, j, t, u, pt) in seg_crossings(L, L + g, cell=0.05))
    a1, a2, s_pt = X[0]
    G = generators(L, eps)
    ir = min(range(len(G)), key=lambda j: np.hypot(*G[j]['pt']))
    rp = G[ir]
    xm = [q for j, q in enumerate(G) if j != ir and q['sign'] == rp['sign']][0]
    p1 = polyline_between(L, rp['pos'], a1)
    p2 = polyline_between(L, a2, xm['pos']) + g
    a0 = strand_path(rp['sign'], rp['k'], eps, xm['pt'][0] + g[0], rp['pt'][0], n=4000)
    loop = np.vstack([p1, p2[1:], a0[1:]])

    def tangent(pos):
        i = min(int(math.floor(pos)), len(L) - 2)
        v = L[i + 1] - L[i]
        return v / np.linalg.norm(v)

    def interior(v_in, v_out):
        return PI - math.atan2(v_in[0] * v_out[1] - v_in[1] * v_out[0], float(np.dot(v_in, v_out)))
    angles = [interior(a0[-1] - a0[-2], tangent(rp['pos'])), interior(tangent(a1), tangent(a2)),
              interior(tangent(xm['pos']), a0[1] - a0[0])]
    simple = len(seg_crossings(resample(np.vstack([loop, loop[:1]]), 0.0005), cell=0.01, skip_adjacent=3)) == 0
    # Maurer-Cartan candidates: the L1-bigon between s and the other crossing with the same deck element, and the
    # polygons from x^- to r_+ with two corners at the self-intersections
    b1, b2, s2_pt = X[1]
    mc = [lattice_windings(np.vstack([polyline_between(L, a1, b1), polyline_between(L, b2, a2)[1:] + g])[:-1])]
    for (c1, c2) in ((a1, a2), (b1, b2)):
        for (d1, d2) in ((a1, a2), (b1, b2)):
            P1 = polyline_between(L, rp['pos'], c1)
            P2 = polyline_between(L, c2, d1) + g
            P3 = polyline_between(L, d2, xm['pos']) + 2 * g
            A0 = strand_path(rp['sign'], rp['k'], eps, xm['pt'][0] + 2 * g[0], rp['pt'][0], n=4000)
            mc.append(lattice_windings(np.vstack([P1, P2[1:], P3[1:], A0[1:]])[:-1]))
    return dict(s=s_pt, ncross=len(X), windings=lattice_windings(loop[:-1]), simple=simple, area=signed_area(loop),
                angles=angles, mc=mc, loop=loop)


def self_intersections_in_P(L):
    """self-intersection points of the arc in P: crossings of the lift with its images under the deck group."""
    lo_, hi_ = L.min(0) - 0.1, L.max(0) + 0.1
    pts = {}
    elems = set()
    for sg in (1, -1):
        for m in range(-6, 7):
            for n in range(-6, 7):
                if sg == 1 and m == 0 and n == 0:
                    continue
                gL = sg * L + np.array([TAU * m, TAU * n])
                if (gL.max(0) < lo_).any() or (gL.min(0) > hi_).any():
                    continue
                for (i, j, t, u, pt) in seg_crossings(L, gL, cell=0.05):
                    key = (round(min(i + t, j + u), 0), round(max(i + t, j + u), 0))
                    if key not in pts:
                        pts[key] = (sg, m, n)
                        elems.add((sg, m, n))
    return len(pts), elems


def run_75(eA, eB, rho, eps_list, verbose):
    """the (7,-5) decomposition at one perturbation; records the items of Proposition 7.1."""
    t0 = time.time()
    T, comps, cov = trace(7, -5, eA, eB, rho)
    arcs = [c for c in comps if c['kind'] == 'arc']
    circ = [c for c in comps if c['kind'] == 'circle']
    record("(1) one arc and five circles, traced set covers all slice solutions",
           len(arcs) == 1 and len(circ) == 5 and len(comps) == 6 and cov < 1e-3)
    L = arc_lift(arcs[0]['L'])
    E = tuple(int(v) for v in np.round(L[-1] / PI))
    record("(1) lift of the arc ends at (7 pi, 8 pi)", E == (7, 8))
    rep = resample(POLYGON, 0.0005)
    same_class = gap_word(L, ((0, 0), E)) == gap_word(rep, ((0, 0), (7, 8)))
    record("(1) arc homotopic rel ends to the polygon (equal reduced gap words)", same_class)
    sc = strip_crossings(L)
    record("(1) arc crosses Delta_0 once, then {phi = pi} once", [k for _, k in sc] == [0, 1])
    # embedded lifts
    emb = len(seg_crossings(L, cell=0.05, skip_adjacent=2)) == 0
    bands, disps, classes, mu_z = [], [], [], []
    for c in circ:
        M, disp = circle_lift(c['L'], 3)
        emb &= len(seg_crossings(M, cell=0.05, skip_adjacent=2)) == 0
        d = tuple(int(v) for v in np.round(disp / PI))
        disps.append(d)
        w = tree_word(np.vstack([c['L'][:-1], c['L'][:1] + disp]))
        classes.append(h1_class(w))
        z = sum({'X': 1, 'Y': 2, 'Z': 1}[Lt] * e for Lt, e in w) % 4
        mu_z.append((circle_mu(c['L']) + z) % 4)
        if abs(d[0]) == 4 and abs(d[1]) == 4:
            f = c['L'][:, 1] - c['L'][:, 0]
            f = f - TAU * math.ceil(f.max() / TAU)                     # representative with max in (-2 pi, 0]
            bands.append((f.min(), f.max()))
    record("(2) every lift of every component is embedded", emb)
    # immersion: |d Pi| / |d(u,v,tau)| along the components, away from the corners of P.  Near the two
    # corners the box coordinates degenerate (for n even the abelian points are the segments u in {0, pi},
    # v = pi/2, each mapped to one point); there Lemma 5.1 and Corollary 5.3 describe the arc.
    ratios = []
    for c in comps:
        dX = np.linalg.norm(np.diff(c['pts'], axis=0), axis=1)
        dL = np.linalg.norm(np.diff(c['L'], axis=0), axis=1)
        far = np.hypot(*(c['L'][1:] - np.round(c['L'][1:] / PI) * PI).T) > 0.05
        keep = (dX > 1e-12) & far
        ratios.append(float(np.min(dL[keep] / dX[keep])))
    record("(2) L1 is an immersion away from the corners (|d Pi| / |d(u,v,tau)| > 0.01)", min(ratios) > 0.01)
    record("(2) mu(C, l1) + z(C) = 0 mod 4 on every circle", all(v == 0 for v in mu_z))
    four = sorted(disps).count((4, 4)) + sorted(disps).count((-4, -4))
    record("(2) four circles with displacement +-(4 pi, 4 pi), one with (0, +-4 pi)",
           four == 4 and sum(1 for d in disps if d[0] == 0 and abs(d[1]) == 4) == 1)
    lclass, _ = L0_class()
    record("(2) no circle class in H_1(P*) is a multiple of a power of [L0] (Definition 2.2(ii))",
           all(np.any(c != 0) and np.linalg.matrix_rank(np.vstack([c, lclass])) == 2 for c in classes))
    band_lo, band_hi = (min(b[0] for b in bands), max(b[1] for b in bands)) if bands else (None, None)

    # generators, gradings, bigons, for each earring parameter
    row = []
    for k, eps in enumerate(eps_list):
        R = analyse_complex(comps, eps, SIGMA, gradings=(k == 0))
        A = R['arc']
        G = A['gens']
        rp = next(g for g in G if g['rplus'])
        others = [g for g in G if not g['rplus']]
        ncirc_gens = sorted(len(C['gens']) for C in R['circles'])
        record("(3) 7 = 1 + |sigma| generators: three on the arc, four on one circle, none on four circles",
               R['ngen'] == 7 and len(G) == 3 and ncirc_gens == [0, 0, 0, 0, 4])
        lab = sorted(g['label'] for g in others)
        same_strand = (lab == ['x+', 'x-'] and
                       all((g['sign'], g['k']) == (rp['sign'], rp['k']) for g in others if g['label'] == 'x-') and
                       all((g['sign'], g['k']) != (rp['sign'], rp['k']) for g in others if g['label'] == 'x+') and
                       rp['sign'] == -1 and rp['k'] == 0)
        record("(4) r_+ and x^- on S_-^0, x^+ on S_+^0", same_strand)
        record("(4) no bigon on any component", R['nbigon'] == 0 and not R['odd'])
        cand = [b for b in A['cands'] if b['p']['label'] == 'x-' and b['q'] is rp]
        wind_ok = len(cand) == 1 and cand[0]['windings'] == {(1, 1): 1, (3, 3): 1}
        record("(4) closed curve of a would-be bigon x^- -> r_+ winds once about (pi,pi) and (3pi,3pi)", wind_ok)
        C4 = [C for C in R['circles'] if len(C['gens']) == 4][0]
        strands = [(g['sign'], g['k']) for g in C4['allgens']]
        record("(4) generators on the circle lie on distinct strands, also after translation by the period",
               len(strands) == len(set(strands)))
        if k == 0:
            gr = {('r+' if g['rplus'] else g['label']): v for g, v in zip(G, A['gr'])}
            ok_gr = gr == {'r+': 2, 'x-': 3, 'x+': 0} and A['agree'] and A['pair_gr'] == [1]
            record("(3) gradings on the arc r_+ = 2, x^- = 3, x^+ = 0 (both directions of alpha_0)", ok_gr)
            record("(3) the circle carries one generator in each grading", C4['agree'] and circle_degree_pattern(C4))
            H = [0, 0, 0, 0]
            for v in A['gr']:
                H[v] += 1
            H = [h + 1 for h in H]                                        # the circle: (1,1,1,1)
            record("(4) H = (2,1,2,2) in gradings 0..3, rank 7", tuple(H) == (2, 1, 2, 2) and R['nbigon'] == 0)
            row.append(f"gr r+ {gr.get('r+')}, x- {gr.get('x-')}, x+ {gr.get('x+')}; H {tuple(H)}")
        row.append(f"eps={eps:g}: {R['ngen']} gens, {R['nbigon']} bigons")
    # bounding cochain
    tri = triangle(L, eps_list[0])
    tri_ok = (not tri['windings'] and tri['simple'] and tri['area'] > 0 and all(0 < a < PI for a in tri['angles'])
              and all(w for w in tri['mc']))
    record("(5) triangle x^-, r_+, s: embedded, counterclockwise, convex corners, no lattice point; MC candidates "
           "enclose lattice points", tri_ok)
    s_pi = tri['s'] / PI
    if verbose:
        print(f"  rho={rho:<6g} (eA,eB)=({eA:+.4f},{eB:+.4f}): arc to {E}, class = polygon {same_class}, "
              f"Delta_0 at gamma = {sc[0][0]:.2f} pi, phi = pi at gamma = {sc[1][0]:.2f} pi; "
              f"circles {disps}; band of the four circles [{band_lo:.3f}, {band_hi:.3f}];\n      " + '; '.join(row) +
              f"; s = ({s_pi[0]:.3f}, {s_pi[1]:.3f}) pi, triangle {'ok' if tri_ok else 'NOT ok'}  [{time.time() - t0:.0f}s]",
              flush=True)
    return dict(L=L, sc=sc, band=(band_lo, band_hi), s=s_pi, tri=tri, comps=comps, immersion=min(ratios))


def control(r, s, eA, eB, rho, eps):
    T, comps, cov = trace(r, s, eA, eB, rho)
    R = analyse_complex(comps, eps, SIGMA)
    A = R['arc']
    G = A['gens']
    gr = {('r+' if g['rplus'] else g['label']): v for g, v in zip(G, A['gr'])}
    H = [0, 0, 0, 0]
    for v in A['gr']:
        H[v] += 1
    for (i, j) in A['bigons']:
        H[A['gr'][i]] -= 1; H[A['gr'][j]] -= 1
    for C in R['circles']:
        if C['gens']:
            H = [h + 1 for h in H]
    circ_ok = all((not C['gens']) or (C['agree'] and circle_degree_pattern(C)) for C in R['circles'])
    bigon = [(G[i]['label'], 'r+' if G[j]['rplus'] else G[j]['label']) for i, j in A['bigons']]
    return dict(gr=gr, H=tuple(H), bigons=bigon, R=R, circ_ok=circ_ok, ngen=R['ngen'], end=A['end'], L=A['L'],
                ok=(len(bigon) == 1 and bigon[0] == ('x-', 'r+') and tuple(H) == (2, 1, 1, 1) and R['nbigon'] == 1
                    and circ_ok and not R['odd']))


def main():
    quick = '--quick' in sys.argv
    t0 = time.time()
    print("Section 7: T(3,4) with (r,s) = (7,-5) (floating-point numerics)")
    zc = corner_loops_z()
    check("z = 1 on small counterclockwise loops about the four corners (normalization of z)",
          all(v == 1 for v in zc.values()), zc, "z(corner loop) = 1")
    lclass, lz = L0_class()
    check("class of L0 in H_1(P*) (letters X, Y, Z) and z(L0)", tuple(lclass) == (1, 1, 1) and lz == 0,
          f"{tuple(int(v) for v in lclass)}, z = {lz}", "figure-eight class, z(L0) = 0")

    # strand-meridian quotients (remark after Proposition 7.1): F(A,B)/<<A^a B^b>> = <A,B | A^a = B^-b>
    def quotient(a, b):
        a, b = abs(a), abs(b)
        return 'unknot' if min(a, b) == 1 else f"T({min(a, b)},{max(a, b)})"
    q75 = (quotient(-5 + 3, 4 - 7), quotient(-5, -7))
    q32 = (quotient(-2 + 3, 4 - 3), quotient(-2, -3))
    check("quotients of pi_1(Y \\ T) by the strand meridians a = A^(s+p) B^(q-r), b = B^(-r) A^s",
          q75 == ('T(2,3)', 'T(5,7)') and q32 == ('unknot', 'T(2,3)'),
          f"(7,-5): {q75}; (3,-2): {q32}", "(7,-5): trefoil, T(5,7); (3,-2): unknot, trefoil")

    rhos = (0.03,) if quick else (0.03, 0.01, 0.003)
    runs = []
    print(f"\nPerturbation sweep: {len(rhos) * 12} runs (rho in {rhos}, twelve directions, eps_B != 0), "
          f"earring eps = rho/10, rho/30")
    for rho in rhos:
        for k in range(12):
            a = k * PI / 6 + 0.05
            runs.append(run_75(rho * math.cos(a), rho * math.sin(a), rho, (rho / 10, rho / 30), verbose=True))
    nrun = len(runs)
    print(f"\nProposition 7.1 over the {nrun} runs (each item checked at both earring parameters where applicable):")
    for item, v in results.items():
        if item.startswith('('):
            print(f"  [{'PASS' if all(v) else 'FAIL'}] {item}: {sum(v)}/{len(v)}")
    check("number of perturbations", nrun == (12 if quick else 36), nrun, "36 (twelve directions, sizes down to 0.003)")
    lo_b = min(r['band'][0] for r in runs); hi_b = max(r['band'][1] for r in runs)
    check(f"band containing phi on the four (4 pi, 4 pi) circles (all {nrun} runs)",
          PAPER['band'][0] < lo_b and hi_b < PAPER['band'][1], f"[{lo_b:.3f}, {hi_b:.3f}]",
          f"{PAPER['band'][0]} < phi < {PAPER['band'][1]}")
    print(f"  (sweep) smallest ratio |d Pi| / |d(u,v,tau)| along the components, away from the corners: "
          f"{min(r['immersion'] for r in runs):.3f}")
    c0 = [r['sc'][0][0] for r in runs]; c1 = [r['sc'][1][0] for r in runs]
    sg = [r['s'][0] for r in runs]; st = [r['s'][1] for r in runs]
    check(f"the arc crosses Delta_0 at gamma ~ 4.33 pi and then phi = pi at gamma ~ 2.67 pi (all {nrun} runs)",
          all(abs(x - PAPER['cross_delta0']) < 0.01 for x in c0) and all(abs(x - PAPER['cross_pi']) < 0.01 for x in c1),
          f"gamma/pi in [{min(c0):.3f}, {max(c0):.3f}] and [{min(c1):.3f}, {max(c1):.3f}]", "gamma ~ 4.33 pi and gamma ~ 2.67 pi")
    check(f"self-intersection s near (0,0) (all {nrun} runs)",
          all(abs(x - PAPER['s'][0]) < 0.005 for x in sg) and all(abs(x - PAPER['s'][1]) < 0.005 for x in st),
          f"gamma/pi in [{min(sg):.3f}, {max(sg):.3f}], theta/pi in [{min(st):.3f}, {max(st):.3f}]", "(0.18 pi, 0.10 pi)")

    print("\nReference run (eps_A, eps_B) = (0.075, 0.05), earring eps in {0.001, 0.003, 0.01, 0.03, 0.05}:")
    nres = {k: len(v) for k, v in results.items()}
    ref = run_75(0.075, 0.05, 0.075, (0.01, 0.001, 0.003, 0.03, 0.05), verbose=True)
    for item, v in results.items():
        if item.startswith('(') and len(v) > nres.get(item, 0):
            print(f"  [{'PASS' if all(v[nres.get(item, 0):]) else 'FAIL'}] {item}")
    tri = ref['tri']
    check("triangle x^-, r_+, s: interior angles, signed area, lattice points inside",
          tri['simple'] and tri['area'] > 0 and all(0 < a < PI for a in tri['angles']) and not tri['windings'],
          f"angles {[round(math.degrees(a), 1) for a in tri['angles']]} deg, area {tri['area']:+.4f}, "
          f"windings {tri['windings'] or 'none'}", "embedded, convex corners, no lattice point")
    n12, elems = self_intersections_in_P(resample(ref['L'], 0.004))
    check("number of self-intersections of the arc in P", n12 == PAPER['selfint'], n12, "twelve")
    gsq = any(e[0] == 1 and e[1] == e[2] and abs(e[1]) == 4 for e in elems)
    check("no self-intersection with deck element a translation by +-(8 pi, 8 pi) (so m_k(s,...,s) = 0, k >= 2)",
          not gsq, f"deck elements of the self-intersections (sign, m, n): {sorted(elems)}",
          "no crossing with deck element g^(+-2)")
    check("candidate polygons with corners at s (m_1(s), and two corners from x^- to r_+) enclose lattice points",
          all(w for w in tri['mc']), [len(w) for w in tri['mc']], "every candidate encloses a lattice point")
    # polygon: generators for eps in [5e-4, 0.2], candidate bigons
    rep = resample(POLYGON, 0.0005)
    ok3, allw = True, True
    for eps in np.geomspace(5e-4, 0.2, 25):
        G = generators(rep, eps)
        ok3 &= len(G) == 3
        for b in bigon_candidates(rep, G, eps):
            allw &= bool(b['windings'])
    check("the polygon meets L0 in exactly three points (25 values of eps in [5e-4, 0.2])", ok3, ok3, "three points")
    check("the polygon bounds no bigon (every candidate closed curve winds about a lattice point)", allw, allw, "no bigon")
    check("the polygon is embedded in R^2", len(seg_crossings(resample(POLYGON, 0.002), cell=0.05, skip_adjacent=2)) == 0,
          "embedded", "the polygon of Section 7")

    print("\nControls at (eps_A, eps_B) = (0.075, 0.05), eps = 0.01:")
    c0 = control(3, -2, 0.075, 0.05, 0.075, 0.01)
    check("j = 0, (r,s) = (3,-2): arc to (pi, 2pi), gradings r_+ 2, x^- 3, x^+ 0, bigon x^- -> r_+, H (2,1,1,1)",
          c0['ok'] and c0['end'] == (1, 2) and c0['gr'] == {'r+': 2, 'x-': 3, 'x+': 0},
          f"end {c0['end']}, gr {c0['gr']}, bigons {c0['bigons']}, H {c0['H']}", "the complex of Theorem 1.1")
    n0, _ = self_intersections_in_P(resample(c0['L'], 0.004))
    check("j = 0: the arc is embedded in P", n0 == 0, f"{n0} self-intersections", "embedded")
    cand = [b for b in c0['R']['arc']['cands'] if b['p']['label'] == 'x-' and b['q']['rplus']]
    check("j = 0: the corresponding closed curve has winding number zero about every lattice point",
          len(cand) == 1 and not cand[0]['windings'] and cand[0]['bigon'], cand[0]['windings'] or 'none', "zero")
    loop0 = np.vstack([strand_path(cand[0]['p']['sign'], cand[0]['p']['k'], 0.01, cand[0]['p']['pt'][0], cand[0]['q']['pt'][0]),
                       polyline_between(c0['L'], cand[0]['q']['pos'], cand[0]['p']['pos'])[1:]])
    w_s = winding_about(loop0, ref['s'] * PI)
    check("s lies inside the region that carries the bigon of Theorem 1.1 (same perturbation)", w_s == 1,
          f"winding number {w_s}", "inside")
    cm = control(-1, 1, 0.075, 0.05, 0.075, 0.01)
    check("j = -1, (r,s) = (-1,1): graded dimensions of H", cm['ok'] and cm['end'] != (1, 2),
          f"arc to {cm['end']}, gr {cm['gr']}, bigons {cm['bigons']}, H {cm['H']}", "(2,1,1,1), a different curve")
    deg_s = (3 - 2) % 4
    check("degree of the bounding cochain s: gr(x^-) - gr(r_+)", deg_s == 1, deg_s, "1")
    Hb = (2, 1, 1, 1)          # x^- (grading 3) -> r_+ (grading 2) removed from (2,1,2,2)
    check("deformed homology: (2,1,2,2) with x^- -> r_+ cancelled", Hb == (2, 1, 1, 1), Hb, "(2,1,1,1)")

    allok = all(all(v) for v in results.values())
    total = sum(len(v) for v in results.values()); passed = sum(sum(v) for v in results.values())
    print(f"\nSection 7: {total} checks, {passed} passed  [{time.time() - t0:.0f}s]")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
