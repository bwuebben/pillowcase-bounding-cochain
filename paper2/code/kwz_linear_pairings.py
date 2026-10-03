"""Linear test curves: Remark 7.7 (the direction of c), Remark 7.8 (curves parallel to c) and Proposition 7.10.
Exact computations over F2.

Remark 7.7.  For each q = 5, 7, ..., 13, all linear curves R (Definition 7.4) of direction (m, n) with |m|, |n| <= 6
or parallel to alpha_q, on the lattice lines of both orbits, with every side pattern of period 2, 4 or 6, together
with the simple closed curves: 1887 curves for each q.  For each q exactly 22 of them have
dim HF(R, BN_q) != dim HF(R, alpha_q) + dim HF(R, (c, J_2)), all parallel to c: the 19 side patterns on the lattice
line through (0,0) and (pi,pi), the two constant patterns on the line through (pi,0) and (0,pi), and the simple
closed curve; none parallel to alpha_q.  The values for c (11 against 13 at q = 7, 7 against 9 at q = 5) and for the
figure-eight around the arc from (0,0) to (pi,pi) (13 against 11, 9 against 7), and the equality for E_{-1}.  The
constant patterns and the simple closed curve are isotopic to c (their encodings are strictly isomorphic to c).
For q = 5, 7 the deformed object (N_q, b + kappa_q) is also compared with BN_q on every test curve.
(Theorem 7.5 itself is proved in the paper; this census supports Remark 7.7.)

Remark 7.8.  c pairs to 11 with BN_7 and to 13 with (N_7, b_S18); (c, J_2) to 22 and 26;
dim HF((c, J_2), (c, J_2)) = 4; the image under Phi of the special component of Kh~(T_7) (the component other than
the figure-eight Kh~(Q_{-r_7})), which is a linear curve of direction c, pairs to 52 with BN_7 and to 44 with
(N_7, b_S18).  At q = 5: 7, 9; 14, 18; 36, 28.

Proposition 7.10.  BN_7 pairs to 9 with E_{-1/2} and to 31 with E_{-3/4}; among the objects allowed by Hypothesis
6.13 (the 45 distinct switches of Proposition 6.10, the nine excluded census smoothings of Proposition 6.11, and the
16 elements of the span of Proposition 6.12) only (N, b_S18) and (N, b_S25) have this pair, and they are strictly
isomorphic; dim HF((c, J_2), (N, b_S18)) = 22 + 4 = 26, whereas dim HF((c, J_2), BN_7) = 22.

Option --quick: Remark 7.7 for q = 5, 7 only.
"""
import argparse, itertools, os, sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def census_remark77(rep, q, check_deformed):
    a = (1 - q, 3 * q - 4)
    wq = (-a[0], -a[1]) if a[0] < 0 else a
    prims = G.primitive_directions(6)
    dirs = list(prims) + ([wq] if wq not in prims else [])
    BN = G.bn_smith(q)
    al = G.alpha(q)
    cJ = G.c_curve(2)
    deformed = []
    if check_deformed:
        Nq, L = O.N(q)
        deformed = [A.deform(Nq, b, O.kappa(q)) for b in O.selected_b(q).values()]
    total = 0
    mism = []
    nobj = 0
    for w in dirs:
        par = 'c' if w[0] + w[1] == 0 else ('alpha_q' if w[0] * a[1] - w[1] * a[0] == 0 else 'generic')
        curves = G.linear_test_curves(w, 4 if w == wq and w not in prims else 6)
        for label, R in curves:
            total += 1
            vb, va, vc = A.hom_dim(R, BN), A.hom_dim(R, al), A.hom_dim(R, cJ)
            if vb != va + vc:
                mism.append((w, par, label, vb, va + vc, R))
            for Y in deformed:
                nobj += A.hom_dim(R, Y) != vb
    rep.check(f'q={q}: number of linear test curves', 1887, total)
    rep.check(f'q={q}: curves pairing differently with BN_q and with alpha_q + (c, J_2)', 22, len(mism))
    rep.check(f'q={q}: ... all parallel to c', True, all(m[1] == 'c' for m in mism))
    rep.check(f'q={q}: ... none parallel to alpha_q', 0, sum(1 for m in mism if m[1] == 'alpha_q'))
    lr = [m for m in mism if m[2][0] == 'L-R']
    nlr = sum(1 for label, R in G.linear_test_curves(G.CDIR if G.CDIR[0] > 0 else (-G.CDIR[0], -G.CDIR[1]), 6)
              if label[0] == 'L-R')
    rep.check(f'q={q}: ... on the line through (0,0) and (pi,pi): all of its side patterns', (19, 19), (len(lr), nlr))
    star = [m for m in mism if m[2][0] == '*-M']
    rep.check(f'q={q}: ... on the line through (pi,0) and (0,pi): the two constant patterns', True,
              len(star) == 2 and all(len(set(m[2][2])) == 1 for m in star))
    sc = [m for m in mism if m[2][0] == 'scc']
    rep.check(f'q={q}: ... and the simple closed curve', 1, len(sc))
    iso = all(A.strictly_isomorphic(m[5], G.c_curve(1)) for m in star + sc) and \
        all(A.strictly_isomorphic(R, G.c_curve(1)) for label, R in
            G.linear_test_curves((1, -1), 2) if label[0] == 'L-R' and label[1] == 2 and len(set(label[2])) == 1)
    rep.check(f'q={q}: the constant patterns on either line and the simple closed curve are isotopic to c', True, iso)
    if check_deformed:
        rep.check(f'q={q}: (N_q, b + kappa_q) and BN_q pair equally with every test curve (mismatches)', 0, nobj,
                  'expected')
    vals = {}
    for m in mism:
        if m[2] == ('scc',):
            vals['c'] = (m[3], m[4])
        if m[2][0] == 'L-R' and m[2][1] == 2 and m[2][2] == (1, -1):
            vals['f8'] = (m[3], m[4])
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args()
    t0 = time.time()
    rep = A.Report('Remarks 7.7, 7.8 and Proposition 7.10')

    rep.section('Remark 7.7')
    exp = {7: {'c': (11, 13), 'f8': (13, 11)}, 5: {'c': (7, 9), 'f8': (9, 7)}}
    for q in ((5, 7) if args.quick else (5, 7, 9, 11, 13)):
        vals = census_remark77(rep, q, q in (5, 7))
        if q in exp:
            rep.check(f'q={q}: c pairs with BN_q and with alpha_q + (c, J_2)', exp[q]['c'], vals.get('c'))
            rep.check(f'q={q}: the figure-eight around the arc (0,0)-(pi,pi) pairs with BN_q and with alpha_q + (c, J_2)',
                      exp[q]['f8'], vals.get('f8'))
            Em1 = G.earring(F(-1))
            rep.check(f'q={q}: E_-1 pairs equally', True,
                      A.hom_dim(Em1, G.bn_smith(q)) == A.hom_dim(Em1, A.union(G.alpha(q), G.c_curve(2))))
    f8LR = [R for label, R in G.linear_test_curves((1, -1), 2) if label == ('L-R', 2, (1, -1))][0]
    rep.check('the figure-eight around the arc (0,0)-(pi,pi) is not an earring: it differs from E_-1, the earring of '
              'its direction', True, not A.strictly_isomorphic(f8LR, G.earring(F(-1))))

    rep.section('Remark 7.8')
    exp8 = {7: (11, 13, 22, 26, 52, 44), 5: (7, 9, 14, 18, 36, 28)}
    cJ = G.c_curve(2)
    rep.check('dim HF((c, J_2), (c, J_2))', 4, A.hom_dim(cJ, cJ))
    for q in (7, 5):
        Nq, L = O.N(q)
        X = A.deform(Nq, list(O.selected_b(q).values())[0])
        BN = G.bn_smith(q)
        comps = O.kht_kh(q)
        f8 = G.kh_rational(-O.r_q(q))
        special = [c for c in comps if not A.strictly_isomorphic(c, f8)]
        rep.check(f'q={q}: Kh~(T_q) has two components, one of them the figure-eight Kh~(Q_-r_q)', (2, 1),
                  (len(comps), len(comps) - len(special)))
        pts = G.linear_loop_points((0, 1), (1, -1), 8, [1, 1, 1, 1, -1, -1, -1, -1])
        Rlin = G.encode(pts, closed=True)
        rep.check(f'q={q}: the special component is a linear curve of direction c (line through (pi,0), (0,pi))',
                  True, len(special) == 1 and A.strictly_isomorphic(special[0], Rlin))
        phiK = G.encode(G.shift(pts), closed=True)
        got = (A.hom_dim(G.c_curve(1), BN), A.hom_dim(G.c_curve(1), X), A.hom_dim(cJ, BN), A.hom_dim(cJ, X),
               A.hom_dim(phiK, BN), A.hom_dim(phiK, X))
        rep.check(f'q={q}: c, (c, J_2) and Phi(special component) against BN_q and (N_q, b)', exp8[q], got)

    rep.section('Proposition 7.10')
    N7, L7 = O.N7()
    E, _ = O.E()
    Es, _ = O.Estar()
    BN7 = G.bn_smith(7)
    rep.check('BN_7 pairs with E_{-1/2} and E_{-3/4}', (9, 31), (A.hom_dim(E, BN7), A.hom_dim(Es, BN7)))
    kh = O.khovanov_ranks()
    rep.check('rk Kh~ of P(-2,3,7) and K_* = K_{-3/4}(7)', (9, 31), (kh[('-1/2', 7)], kh[('-3/4', 7)]))
    allowed = {}
    census = []
    for line in open(os.path.join(O.DATA, 'q7_census.txt'), encoding='utf-8'):
        if line.startswith('#') or not line.strip():
            continue
        name, kind, gp, sw = line.rstrip('\n').split('\t')
        census.append((name, gp == 'yes', sw))
    distinct = {}
    for name, gp, sw in census:
        if gp:
            distinct.setdefault(frozenset(O._terms(sw, O.N7_LABELS)), []).append(name)
    for b, names in distinct.items():
        allowed['switch ' + '/'.join(names)] = A.deform(N7, list(b))
    smooth = O.complexes('q7_smoothings')
    excluded = [name for name, gp, sw in census if not gp]
    for name in excluded:
        allowed['excluded ' + name] = smooth[name + '_census'][0]
    named = ['S18', 'S25', 'S69', 'S74']
    for mask in range(16):
        sel = [named[i] for i in range(4) if mask >> i & 1]
        b = []
        for s in sel:
            b = A.xor(b, O.switch7(s))
        allowed['span ' + ('+'.join(sel) or '0')] = A.deform(N7, b)
    rep.check('allowed objects: distinct switches, excluded orbits, span elements', (45, 9, 16),
              (len(distinct), len(excluded), 16))
    pair931 = sorted(k for k, X in allowed.items() if A.is_mc(X) and A.hom_dim(E, X) == 9 and A.hom_dim(Es, X) == 31)
    rep.check('allowed objects with the pair (9, 31)', ['span S18', 'span S25', 'switch S18', 'switch S25'], pair931)
    XS18 = A.deform(N7, O.switch7('S18'))
    rep.check('(N, b_S18) and (N, b_S25) are strictly isomorphic', True,
              A.strictly_isomorphic(XS18, A.deform(N7, O.switch7('S25'))))
    al7 = G.alpha(7)
    rep.check('dim HF((c, J_2), alpha_7), dim HF((c, J_2), (c, J_2)), dim HF((c, J_2), (N, b_S18))', (22, 4, 26),
              (A.hom_dim(cJ, al7), A.hom_dim(cJ, cJ), A.hom_dim(cJ, XS18)))
    rep.check('dim HF((c, J_2), BN_7)', 22, A.hom_dim(cJ, BN7))
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
