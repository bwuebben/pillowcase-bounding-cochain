"""Theorem 7.1 (the structure theorem), Remark 7.3 and Appendix A.1.

Exact computations in Tw(B) over F2.  Every statement of Theorem 7.1 asserts that two explicit reduced twisted
complexes are strictly isomorphic, or that an explicit matrix over B squares to zero; all complexes involved are of
loop type, so strict isomorphism is decided by comparing the words of their components.

  (a) (N_q, b) is strictly isomorphic to alpha_q + (c, J_2), for q = 7, b in {b_S18, b_S25}, and q = 5, b = b_5.
  (b) For odd 5 <= q <= 21, BN_q = Phi(BN~(T_q)) is the arc obtained from alpha_q by sliding its end at (0,0) twice
      around that corner, parallel to c; among the 24 slides of alpha_q it is the only one.
  (c) (delta + kappa_q)^2 = (delta + b + kappa_q)^2 = 0 and (N_q, b + kappa_q) is strictly isomorphic to BN_q.
Remark 7.3: for t > 0 the comparison object is BN~(T_q) itself, the slide is at (pi,0), and (c) holds with kappa'_q;
(N_q, b + kappa'_q) pairs to 13 (q = 7) and 9 (q = 5) with E_{-1}, against rk Kh~(K_{-1}(q)) = 11 and 7.
Appendix A.1: the complex N_5, the smoothing b_5 and the corner term kappa_5.
Section 7.2: alpha_q and BN_q are bigraded arcs; c and (c, J_2) admit a delta-grading but no bigrading.
BN~(T_q) is read from the output of kht++ in data/kht (see README for how to regenerate it).
"""
import sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def end_at_origin(X, q):
    """For an object X whose arc component is strictly isomorphic to alpha_q: {arc end of X: '(0,0)' or '(pi,0)'}.
    The encoding of alpha_q = [R0, M_q] starts at its end over (0,0)."""
    al = G.alpha(q)
    (seqA, labsA, _), = A.walk(al)
    assert seqA[0] == 0 and G.corner_of(G.R0) == '(0,0)'
    wordA = [al[0][seqA[0]]] + [x for t in range(len(labsA)) for x in (labsA[t], al[0][seqA[t + 1]])]
    out = {}
    for seq, labs, closed in A.walk(X):
        if closed:
            continue
        for orient in (1, -1):
            s_ = seq if orient == 1 else seq[::-1]
            l_ = labs if orient == 1 else [A._flip(l) for l in reversed(labs)]
            w = [X[0][s_[0]]] + [x for t in range(len(l_)) for x in (l_[t], X[0][s_[t + 1]])]
            if w == wordA:
                out[s_[0]] = '(0,0)'
                out[s_[-1]] = '(pi,0)'
    return out


def grading_obstructions(X):
    """For each closed component of a loop-type object, the obstruction (h, q, delta) to a bigrading: the sum around
    the loop of (1, 0) - gr(a) over its arrows a (with sign by direction), with the gradings gr D = (0, -2),
    gr S = (0, -1) of Kotelskiy, Watson and Zibrowius and delta = q/2 - h.  A component is bigradable iff (h, q) = 0
    and delta-gradable iff delta = 0."""
    out = []
    for seq, labs, closed in A.walk(X):
        if not closed:
            continue
        h = qq = 0
        for kd, k, dr in labs:
            qa = -2 * k if kd == 'D' else -k
            h += dr
            qq -= dr * qa
        out.append((h, qq, F(qq, 2) - h))
    return out


def main():
    t0 = time.time()
    rep = A.Report('Theorem 7.1, Remark 7.3 and Appendix A.1')
    pl = O.complexes('pl_models')
    traced = O.complexes('traced_curves')
    E, _ = O.E()
    Es, _ = O.Estar()
    E_half, E_34, E_m1 = G.earring(F(-1, 2)), G.earring(F(-3, 4)), G.earring(F(-1))
    kh = O.khovanov_ranks()

    rep.section('the curves of Section 7.2 and the translation Phi')
    for q in (5, 7):
        a = G.alpha(q)
        rep.check(f'q={q}: Phi(alpha_q) strictly isomorphic to alpha_q', True,
                  A.strictly_isomorphic(G.encode(G.shift(G.alpha_points(q))), a))
    rep.check('Phi(c) strictly isomorphic to c', True,
              A.strictly_isomorphic(G.encode(G.shift(G.c_points(1)), closed=True), G.c_curve(1)))
    rep.check('Phi(doubled c) strictly isomorphic to the doubled word of c', True,
              A.strictly_isomorphic(G.encode(G.shift(G.c_points(2)), closed=True), G.c_curve(2)))
    rep.check('alpha_q joins (0,0) and (pi,0) (q = 5, 7)', True,
              all({G.corner_of(G.R0), G.corner_of(G.Mq(q))} == {'(0,0)', '(pi,0)'} for q in (5, 7)))

    rep.section('gradings (Section 7.2)')
    ob_c, ob_cJ = grading_obstructions(G.c_curve(1)), grading_obstructions(G.c_curve(2))
    rep.check('c admits a delta-grading but no bigrading', True, len(ob_c) == 1 and ob_c[0][2] == 0 and ob_c[0][:2] != (0, 0))
    rep.check('(c, J_2) (the doubled word) admits a delta-grading but no bigrading', True,
              len(ob_cJ) == 1 and ob_cJ[0][2] == 0 and ob_cJ[0][:2] != (0, 0))
    rep.note(f'obstructions (h, q, delta): c {ob_c}, doubled c {ob_cJ}')
    rep.check('alpha_q and BN_q are arcs (no obstruction), q = 5, 7, ..., 21', True,
              all(not grading_obstructions(G.alpha(q)) and not grading_obstructions(G.bn_smith(q)) for q in range(5, 23, 2)))
    rs = (F(-1, 2), F(-3, 4), F(3, 4), F(-1), F(1, 3), O.r_q(5), O.r_q(7))
    rep.check('the earrings are bigraded after moving back by Phi: Kh~(Q_-r) has no obstruction (r = -1/2, -3/4, 3/4, '
              '-1, 1/3, r_5, r_7); E_r itself is delta-graded', True,
              all(o[:2] == (0, 0) for r in rs for o in grading_obstructions(G.kh_rational(-r)))
              and all(o[2] == 0 for r in rs for o in grading_obstructions(G.earring(r))), 'expected')

    rep.section('Appendix A.1: the complex N_5')
    N5, L5 = O.N5()
    rep.check('N_5: number of generators', 23, len(N5[0]))
    rep.check('N_5: iota_circ generators', [3, 4, 9, 10, 13, 14, 19, 20], [L5[i] for i in range(23) if N5[0][i] == 'w'])
    rep.check('N_5: delta^2 = 0', True, A.is_mc(N5))
    starts = {}
    for (i, kd, k, j) in N5[1]:
        starts.setdefault(i, set()).add(A.label_of(N5, (i, kd, k, j))[0])
    rep.check('N_5: no arrow ends where an arrow with the same face label begins', True,
              all(A.label_of(N5, t)[0] not in starts.get(t[3], set()) for t in N5[1]))
    comps = A.components(N5)
    rep.check('N_5 is a single arc: number of components', 1, len(comps))
    ends = end_at_origin(A.deform(N5, O.b5()), 5)
    rep.check('N_5: ends (label at (0,0), label at (pi,0))', (0, 22),
              (L5[[e for e, c in ends.items() if c == '(0,0)'][0]], L5[[e for e, c in ends.items() if c == '(pi,0)'][0]]))
    rep.check('N_5 strictly isomorphic to the encoding of L_{5,PL}', True, A.strictly_isomorphic(N5, pl['N5'][0]))
    t5 = [k for k in traced if k.startswith('q5_')]
    rep.check('the traced curves at q = 5 strictly isomorphic to N_5 (number of traces)', len(t5),
              sum(A.strictly_isomorphic(traced[k][0], N5) for k in t5), 'expected')
    rep.check('N_5: pairings with E_{-1/2} and E_{-3/4}', (5, 15), (A.hom_dim(E_half, N5), A.hom_dim(E_34, N5)))
    b5 = O.b5()
    X5 = A.deform(N5, b5)
    rep.check('b_5 exchanges the targets of 8 -R-> 7 and 16 -R-> 15', True,
              sorted(A.arrow_str(N5, t, L5) for t in b5) == sorted(['8 -R-> 7', '16 -R-> 15', '8 -R-> 15', '16 -R-> 7']))
    rep.check('(N_5, b_5): D(b_5) = 0 and b_5^2 = 0, i.e. Maurer--Cartan', True, A.is_mc(X5) and A.is_mc((N5[0], b5)))
    cs = A.components(X5)
    sets = []
    for seq, labs, closed in A.walk(X5):
        sets.append((closed, sorted(L5[v] for v in seq)))
    arcset = [s for c, s in sets if not c]
    loopset = [s for c, s in sets if c]
    rep.check('(N_5, b_5): the arc', [list(range(0, 8)) + list(range(16, 23))], arcset)
    rep.check('(N_5, b_5): the closed component', [list(range(8, 16))], loopset)
    loop = [c for c in cs if len(c[0]) == 8][0]
    rep.check('(N_5, b_5): the closed component is the doubled word of c', True, A.strictly_isomorphic(loop, G.c_curve(2)))
    rep.check('(N_5, b_5) strictly isomorphic to the census smoothings of L_{5,PL} at two orbits', (True, True),
              (A.strictly_isomorphic(X5, pl['N5_orb13'][0]), A.strictly_isomorphic(X5, pl['N5_orb18'][0])))
    rep.check('(N_5, b_5): pairing with E_{-1/2}', 7, A.hom_dim(E_half, X5))
    Y5 = A.deform(N5, b5, O.kappa(5))
    (seq, labs, closed), = A.walk(Y5)
    order = [L5[v] for v in seq]
    target = [22, 21, 20, 19, 18, 17, 16, 7, 6, 5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 15, 14, 13, 12]
    rep.check('(N_5, b_5 + kappa_5) is the single arc 22,...,16,7,...,0,11,10,9,8,15,...,12', True,
              not closed and (order == target or order == target[::-1]))
    rep.check('(N_5, b_5 + kappa_5): pairings with E_{-1/2} and E_{-3/4}', (7, 21),
              (A.hom_dim(E_half, Y5), A.hom_dim(E_34, Y5)))
    rep.check('rk Kh~ of P(-2,3,5) and of K_{-3/4}(5)', (7, 21), (kh[('-1/2', 5)], kh[('-3/4', 5)]))

    rep.section('Section 6.3: N_7 and the encodings of L_{7,PL}')
    N7, L7 = O.N7()
    rep.check('N_7 strictly isomorphic to the encoding of L_{7,PL}', True, A.strictly_isomorphic(N7, pl['N7'][0]))
    rep.check('(N_7, b_S18), (N_7, b_S25) strictly isomorphic to the encoded census smoothings at S_18, S_25',
              (True, True), (A.strictly_isomorphic(A.deform(N7, O.switch7('S18')), pl['N7_orb18'][0]),
                             A.strictly_isomorphic(A.deform(N7, O.switch7('S25')), pl['N7_orb25'][0])))
    ends7 = end_at_origin(A.deform(N7, O.switch7('S18')), 7)
    rep.check('N_7: ends (label at (0,0), label at (pi,0))', (36, 18),
              (L7[[e for e, c in ends7.items() if c == '(0,0)'][0]], L7[[e for e, c in ends7.items() if c == '(pi,0)'][0]]))

    rep.section('Theorem 7.1(a)')
    for q in (7, 5):
        Nq, L = O.N(q)
        for bname, b in O.selected_b(q).items():
            X = A.deform(Nq, b)
            sizes = tuple(sorted(A.component_sizes(X), reverse=True))
            rep.check(f'q={q}, {bname}: Maurer--Cartan', True, A.is_mc(X))
            rep.check(f'q={q}, {bname}: component sizes', (23, 8) if q == 7 else (15, 8), sizes)
            arc_c = [c for c in A.components(X) if len(c[0]) == sizes[0]][0]
            loop_c = [c for c in A.components(X) if len(c[0]) == 8][0]
            rep.check(f'q={q}, {bname}: the arc is alpha_q', True, A.strictly_isomorphic(arc_c, G.alpha(q)))
            rep.check(f'q={q}, {bname}: the closed component is the doubled word of c', True,
                      A.strictly_isomorphic(loop_c, G.c_curve(2)))
            rep.check(f'q={q}, {bname}: (N_q, b) strictly isomorphic to alpha_q + (c, J_2)', True,
                      A.strictly_isomorphic(X, A.union(G.alpha(q), G.c_curve(2))))

    rep.section('Theorem 7.1(b), q = 5, 7, ..., 21')
    nq = 0
    for q in range(5, 23, 2):
        BNk = O.kht_bn(q)
        poly_ok = A.strictly_isomorphic(G.encode(G.bn_kwz_points(q)), BNk)
        BN = G.bn_smith(q)
        sl = G.slides(q)
        hits = [k for k, X in sl.items() if A.strictly_isomorphic(X, BN)]
        ok = poly_ok and len(sl) == 24 and hits == [('(0,0)', '-c', 2, 'right')]
        nq += ok
        rep.check(f'q={q}: BN~(T_q) from kht++ is alpha_q slid at (pi,0) in the KWZ chart; the slides of alpha_q that '
                  f'equal BN_q (of 24)', [('(0,0)', '-c', 2, 'right')], hits if poly_ok and len(sl) == 24 else None)
    rep.check('values of q for which (b) holds', 9, nq)
    # geometric meaning of the matching slide: direction -c from p_0 makes an acute angle with p_0 - p_M, and the
    # side 'right' of that direction is the half-plane l^+ not containing p_M
    d = (-G.CDIR[0], -G.CDIR[1])
    acute = all((G.R0[0] - G.Mq(q)[0]) * d[0] + (G.R0[1] - G.Mq(q)[1]) * d[1] > 0 for q in range(5, 23, 2))
    pm_left = all(((G.Mq(q)[0] - G.R0[0]) * d[1] - (G.Mq(q)[1] - G.R0[1]) * d[0]) < 0 for q in range(5, 23, 2))
    rep.check('the slide runs along p_1 = p_0 + (1,-1) (acute angle with p_0 - p_M) on the side l^+ away from p_M',
              True, acute and pm_left)
    rep.check('two periods of c end at p_4 = p_0 + 2(p_2 - p_0), a lift of (0,0)', '(0,0)',
              G.corner_of((G.R0[0] + 4 * d[0], G.R0[1] + 4 * d[1])))

    rep.section('Theorem 7.1(c)')
    for q in (7, 5):
        Nq, L = O.N(q)
        kap = O.kappa(q)
        el = O.end_labels(q)['(0,0)']
        deleted = [t for t in kap if t in Nq[1]]
        added = [t for t in kap if t not in Nq[1]]
        rep.check(f'q={q}: kappa_q deletes one arrow of delta_N (an R-chord) and adds an R-arrow at the end over (0,0)',
                  True, len(deleted) == 1 and len(added) == 1 and A.label_of(Nq, deleted[0]) == 'R'
                  and A.label_of(Nq, added[0]) == 'R' and el in (L[added[0][0]], L[added[0][3]]))
        rep.check(f'q={q}: (delta_N + kappa_q)^2 = 0', True, A.is_mc(A.deform(Nq, kap)))
        for bname, b in O.selected_b(q).items():
            Y = A.deform(Nq, b, kap)
            rep.check(f'q={q}, {bname}: (delta_N + b + kappa_q)^2 = 0', True, A.is_mc(Y))
            rep.check(f'q={q}, {bname}: (N_q, b + kappa_q) strictly isomorphic to BN_q', True,
                      A.strictly_isomorphic(Y, G.bn_smith(q)))

    rep.section('Remark 7.3 (t > 0)')
    for q in (7, 5):
        Nq, L = O.N(q)
        BNk = O.kht_bn(q)
        hits = [k for k, X in G.slides(q).items() if A.strictly_isomorphic(X, BNk)]
        rep.check(f'q={q}: BN~(T_q) = Phi(BN_q) is the slide of alpha_q at (pi,0) (unique among 24)',
                  [('(pi,0)', '+c', 2, 'right')], hits)
        kp = O.kappa_prime(q)
        added = [t for t in kp if t not in Nq[1]]
        rep.check(f"q={q}: kappa'_q adds an arrow at the end over (pi,0) (label {O.end_labels(q)['(pi,0)']})", True,
                  len(added) == 1 and O.end_labels(q)['(pi,0)'] in (L[added[0][0]], L[added[0][3]]))
        for bname, b in O.selected_b(q).items():
            Y = A.deform(Nq, b, kp)
            rep.check(f"q={q}, {bname}: (delta + kappa'_q)^2 = (delta + b + kappa'_q)^2 = 0", True,
                      A.is_mc(A.deform(Nq, kp)) and A.is_mc(Y))
            rep.check(f"q={q}, {bname}: (N_q, b + kappa'_q) strictly isomorphic to BN~(T_q)", True,
                      A.strictly_isomorphic(Y, BNk))
            rep.check(f"q={q}, {bname}: (N_q, b + kappa'_q) is not strictly isomorphic to BN_q", True,
                      not A.strictly_isomorphic(Y, G.bn_smith(q)))
        b = list(O.selected_b(q).values())[0]
        Yp = A.deform(Nq, b, kp)
        rep.check(f"q={q}: dim HF(E_-1, (N_q, b + kappa'_q)) and rk Kh~(K_-1(q))", (13, 11) if q == 7 else (9, 7),
                  (A.hom_dim(E_m1, Yp), kh[('-1', q)]))
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
