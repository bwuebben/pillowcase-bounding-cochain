#!/usr/bin/env python3
r"""Remark 6.12 (and Remark 5.8): numerical computation of the perturbed curves and of the complex C(L0, L1)
for T(3, n), n = 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 22, 23, in the decomposition of Theorem 1.1, at

    (n eps_A, eps_B, eps) = (0.3, 0.05, 0.01)  and  (-0.3, 0.03, 0.02),

with the grading formula (2.2) evaluated directly for every pair of generators on a component.  For each of
the 24 cases the following are compared with the statements of Section 6:

  * the arc: its lift ends at (pi, 2 pi j) and crosses the lines Delta_k listed in Table 1 (Lemma 6.1); it is
    embedded and its interior lies in the open strip 0 < gamma < pi;
  * which resolution occurs (Remark 5.8: split for n = 4, joined for n = 5 (mod 6));
  * the generators: 1 + 2 c_A on the arc, four on every circle, 1 + |sigma| in total (Proposition 6.3);
  * the bigons: none for n = 1, 2 (mod 6), exactly one, from x_1^- to r_+, for n = 4, 5 (Proposition 6.5);
  * the degrees deg = gr - sigma of the arc generators (Proposition 6.9), with gr(x^+, x^-) = 1 at every crossing;
    on every circle one generator in each degree, the two x^- differing by 2;
  * the graded homology against Table 2 and (2.5);
  * for n = 4, 5, 7, the complexes of [HHK18, Sec. 11.6, 11.7, 11.2];
  * for the first perturbation, that the tangent line of L1 is parallel to the line field of slope one only at
    the critical points of phi (two on the arc for n = 4, 5 (mod 6), none otherwise, none on the circles).

Floating-point numerics; these computations are not part of the proofs.
"""
import sys
import time

import numpy as np

from t3n_pillowcase import (PI, TAU, Tangle, decomposition, signature, alexander_norm, inat_degrees,
                            trace_variety, analyse_complex, circle_degree_pattern, seg_crossings,
                            direction_changes, tangent_passages)

NS = (4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 22, 23)
PERTURBATIONS = ((0.3, 0.05, 0.01), (-0.3, 0.03, 0.02))          # (n eps_A, eps_B, eps)

# Table 1 and Proposition 6.9: (n mod 6, resolution) -> (j, lines crossed, degrees (x_i^-, x_i^+) along the arc)
TABLE = {
    (1, None): (1, [], []),
    (2, None): (2, [1], [(1, 2)]),
    (4, 'split'): (1, [0], [(1, 2)]),
    (4, 'joined'): (3, [0, 1, 2], [(1, 2), (3, 0), (1, 2)]),
    (5, 'split'): (2, [0, 1], [(1, 2), (3, 0)]),
    (5, 'joined'): (4, [0, 1, 2, 3], [(1, 2), (3, 0), (1, 2), (3, 0)]),
}
EXPECTED_RESOLUTION = {4: 'split', 5: 'joined'}                    # observed pattern stated in Remark 5.8
HHK_SEC11 = {4: (7, 1, 5), 5: (9, 1, 7), 7: (9, 0, 9)}            # generators, bigons, rank of homology

results = []


def run(n, cA, eB, eps, tangency):
    t0 = time.time()
    r, s = decomposition(n)
    sig = signature(3, n)
    T = Tangle(3, n, r, s, cA / n, eB)
    comps, _ = trace_variety(T)
    R = analyse_complex(comps, eps, sig)
    A = R['arc']
    G = A['gens']
    cls = n % 6
    fails = []

    # arc: crossings with the lines Delta_k, endpoint, embeddedness, open strip
    xs = [g for g in G if not g['rplus']]
    levels = sorted(set(int(round((g['pt'][1] - g['pt'][0]) / TAU)) for g in xs))
    resolution = None
    if cls == 4:
        resolution = 'split' if len(G) == 3 else ('joined' if len(G) == 7 else '?')
    if cls == 5:
        resolution = 'split' if len(G) == 5 else ('joined' if len(G) == 9 else '?')
    j, lines, degs = TABLE.get((cls, resolution), (None, None, None))
    if j is None:
        fails.append(f"arc with {len(G)} generators matches no row of Table 1")
    else:
        if A['end'] != (1, 2 * j):
            fails.append(f"arc ends at {A['end']} (pi units), Table 1: (1, {2 * j})")
        if levels != lines:
            fails.append(f"lines crossed {levels}, Table 1: {lines}")
    if cls in EXPECTED_RESOLUTION and resolution != EXPECTED_RESOLUTION[cls]:
        fails.append(f"resolution {resolution}, Remark 5.8: {EXPECTED_RESOLUTION[cls]}")
    L = A['L']
    inner = L[len(L) // 100: -len(L) // 100, 0]
    if not (np.all(inner > 0) and np.all(inner < PI)):
        fails.append("arc leaves the open strip 0 < gamma < pi")
    if seg_crossings(L, cell=0.05, skip_adjacent=2):
        fails.append("lift of the arc not embedded")

    # generators
    ncirc = len(R['circles'])
    if any(len(C['gens']) != 4 for C in R['circles']):
        fails.append(f"circle generators {[len(C['gens']) for C in R['circles']]}")
    if R['ngen'] != 1 + abs(sig):
        fails.append(f"{R['ngen']} generators, 1 + |sigma| = {1 + abs(sig)}")

    # bigons
    rp = next(i for i, g in enumerate(G) if g['rplus'])
    want_bigons = [] if cls in (1, 2) else 'x1- -> r+'
    if cls in (1, 2):
        if R['nbigon']:
            fails.append(f"{R['nbigon']} bigons, Proposition 6.5: none")
    else:
        ok = (R['nbigon'] == 1 and len(A['bigons']) == 1 and A['bigons'][0][1] == rp
              and any(G[A['bigons'][0][0]] is g for g in xs[:2]) and G[A['bigons'][0][0]]['label'] == 'x-')
        if not ok:
            fails.append(f"bigons {A['bigons']} (arc), {R['nbigon']} in total; Proposition 6.5: exactly x_1^- -> r_+")
    if R['odd']:
        fails.append(f"{len(R['odd'])} null-homotopic candidate loops that are not simple")

    # gradings
    deg = [(gr - sig) % 4 for gr in A['gr']]
    if not A['agree']:
        fails.append("the two directions of alpha_0 give different gradings")
    if any(v != 1 for v in A['pair_gr']):
        fails.append(f"gr(x^+, x^-) = {A['pair_gr']}")
    if deg[rp] != 0:
        fails.append(f"deg r_+ = {deg[rp]}")
    got_pairs = []
    for k in range(0, len(xs) - 1, 2):
        a, b = xs[k], xs[k + 1]
        xm, xp = (a, b) if a['label'] == 'x-' else (b, a)
        got_pairs.append((deg[G.index(xm)], deg[G.index(xp)]))
    if degs is not None and got_pairs != degs:
        fails.append(f"arc degrees {got_pairs}, Proposition 6.9: {degs}")
    if not all(C['agree'] and circle_degree_pattern(C) for C in R['circles']):
        fails.append("circle gradings")
    H = [0, 0, 0, 0]
    for d in deg:
        H[d] += 1
    if cls in (4, 5):
        H[1] -= 1; H[0] -= 1                           # the bigon x_1^- -> r_+ cancels degrees 1 and 0
    for _ in R['circles']:
        H = [h + 1 for h in H]
    if tuple(H) != inat_degrees(3, n):
        fails.append(f"H in degrees 0..3 = {tuple(H)}, (2.5): {inat_degrees(3, n)}")
    if sum(H) != alexander_norm(3, n):
        fails.append("rank of H")
    if n in HHK_SEC11 and (R['ngen'], R['nbigon'], sum(H)) != HHK_SEC11[n]:
        fails.append(f"(generators, bigons, rank) = {(R['ngen'], R['nbigon'], sum(H))}, [HHK18, Sec. 11]: {HHK_SEC11[n]}")

    # tangency of L1 with the line field of slope one (first perturbation)
    tang = ''
    if tangency:
        phi = L[:, 1] - L[:, 0]
        turns, passes = direction_changes(phi), tangent_passages(L)
        circ = [tangent_passages(np.vstack([C['M']])) for C in R['circles']]
        want = 2 if cls in (4, 5) else 0
        if not (turns == passes == want and all(c == 0 for c in circ)):
            fails.append(f"tangency: arc passages {passes}, critical points of phi {turns}, circles {circ}")
        tang = f"; slope-one tangencies on the arc {passes} (= critical points of phi {turns}), on circles {sum(circ)}"

    desc = (f"T(3,{n:2d}) (r,s)=({r},{s}) (n eA, eB, eps)=({cA:+.1f},{eB},{eps}): arc to (pi,{A['end'][1]}pi), "
            f"{resolution or '-'}, crosses Delta_{levels}; {R['ngen']} generators (1+|sigma|={1 + abs(sig)}), "
            f"{ncirc} circles; bigons {'x1- -> r+' if R['nbigon'] == 1 else R['nbigon']}; "
            f"deg r+ {deg[rp]}, (x-,x+) {got_pairs}; H {tuple(H)}{tang}  [{time.time() - t0:.0f}s]")
    ok = not fails
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)
    for f in fails:
        print(f"        mismatch: {f}")


def main():
    t0 = time.time()
    ns = [int(a) for a in sys.argv[1:]] or list(NS)
    print("Remark 6.12: perturbed complexes of T(3,n) (floating-point numerics)")
    for n in ns:
        for k, (cA, eB, eps) in enumerate(PERTURBATIONS):
            run(n, cA, eB, eps, tangency=(k == 0))
    print("\nRemark 6.12: %d cases, %d agree with Theorem 6.10, Table 2 and Remark 5.8  [%.0fs]"
          % (len(results), sum(results), time.time() - t0))
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
