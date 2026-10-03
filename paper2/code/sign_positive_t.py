"""Remark 6.17: the pairings for small t > 0.  Exact computations over F2.

For t > 0 Smith's curve, its self-intersections and their smoothings are the translates by Phi of those for t < 0
(Remark 6.2), while the earring curves do not depend on t.  Since Phi is a symplectomorphism of order two, the
pairings for t > 0 are those of the present objects with the translated earrings Phi(E_r) (which equal KWZ's
Kh~(Q_{-r}), Section 7.1).  Remark 6.17 states:
  * the undeformed pairings remain 7 and 25 (first and second closure);
  * the smoothings at S_18 and S_25 still give (9, 31), those at S_69 and S_74 give (5, 25) and (7, 25);
  * dimension nine against the first earring is reached by seven of the 45 distinct switches, those at S_10, S_18,
    S_22, S_25, S_29, S_30 and the one shared by S_26, S_27, S_28, S_32, S_33;
  * among the 82 census smoothings, the 82 opposite smoothings and the 16 elements of the span of Proposition 6.12,
    the pair (9, 31) occurs only at S_18 and S_25; the opposite smoothings there are separated by the third closure
    (69 against 103); in the span only b_S18 and b_S25 reach dimension nine;
  * so Proposition 6.7 does not hold for t > 0: only two of the four objects have dimension nine.
The census smoothings and the opposite smoothings are read from data/q7_smoothings.txt.  The same quantities for
t < 0 (Propositions 6.7 and 6.10 to 6.12) are printed for comparison.
"""
import os, sys, time
from collections import Counter
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def main():
    t0 = time.time()
    rep = A.Report('Remark 6.17: the pairings for t > 0')
    N7, L7 = O.N7()
    slopes = (F(-1, 2), F(-3, 4), F(3, 4))
    neg = {r: G.earring(r) for r in slopes}                                         # t < 0
    pos = {r: G.encode(G.shift(G.earring_points(r)), closed=True) for r in slopes}  # t > 0: Phi(E_r)
    E, _ = O.E()
    Es, _ = O.Estar()
    rep.check('E_{-1/2}, E_{-3/4} are the complexes (6.2), (6.3)', (True, True),
              (A.strictly_isomorphic(neg[F(-1, 2)], E), A.strictly_isomorphic(neg[F(-3, 4)], Es)))
    rep.check('Phi(E_r) = Kh~(Q_-r) for r = -1/2, -3/4, 3/4', True,
              all(A.strictly_isomorphic(pos[r], G.kh_rational(-r)) for r in slopes))

    def pr(X, sign):
        Ee = neg if sign < 0 else pos
        return (A.hom_dim(Ee[F(-1, 2)], X), A.hom_dim(Ee[F(-3, 4)], X))

    rep.section('the undeformed object and the four named switches')
    rep.check('N: (first, second closure) for t > 0', (7, 25), pr(N7, +1))
    rep.check('N: (first, second closure) for t < 0', (7, 25), pr(N7, -1), 'expected')
    exp_pos = {'S18': (9, 31), 'S25': (9, 31), 'S69': (5, 25), 'S74': (7, 25)}
    exp_neg = {'S18': (9, 31), 'S25': (9, 31), 'S69': (9, 23), 'S74': (9, 25)}
    for s in ('S18', 'S25', 'S69', 'S74'):
        X = A.deform(N7, O.switch7(s))
        rep.check(f'(N, b_{s}) for t > 0', exp_pos[s], pr(X, +1))
        rep.check(f'(N, b_{s}) for t < 0 (Propositions 6.7, 6.8)', exp_neg[s], pr(X, -1), 'expected')
    rep.check('objects of Proposition 6.7 with dimension nine for t > 0', 2,
              sum(pr(A.deform(N7, O.switch7(s)), +1)[0] == 9 for s in exp_pos))

    rep.section('the 45 distinct switches')
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
    rep.check('distinct switches', 45, len(distinct))
    nine_pos = []
    vals = {}
    for b, names in distinct.items():
        X = A.deform(N7, list(b))
        vals['/'.join(names)] = (pr(X, +1), pr(X, -1))
        if vals['/'.join(names)][0][0] == 9:
            nine_pos.append('/'.join(names))
    key = lambda n: int(n.split('/')[0][1:])
    rep.check('switches reaching dimension nine against the first earring for t > 0',
              ['S10', 'S18', 'S22', 'S25', 'S26/S27/S28/S32/S33', 'S29', 'S30'], sorted(nine_pos, key=key))

    rep.section('the 82 census smoothings, the 82 opposite smoothings, the span of Proposition 6.12')
    smooth = O.complexes('q7_smoothings')
    sw_of = {name: frozenset(O._terms(sw, O.N7_LABELS)) for name, gp, sw in census if gp}
    rep.check('the encoded census smoothing equals (N, b_s) at every generator-preserving orbit (number of orbits)', 73,
              sum(A.strictly_isomorphic(smooth[n + '_census'][0], A.deform(N7, list(b))) for n, b in sw_of.items()),
              'expected')
    hits = {+1: [], -1: []}
    third = {}
    for name, gp, sw in census:
        for tag in ('census', 'opposite'):
            X = smooth[f'{name}_{tag}'][0]
            for sign in (+1, -1):
                if pr(X, sign) == (9, 31):
                    hits[sign].append(f'{name} {tag}')
            if name in ('S18', 'S25'):
                third[(name, tag)] = A.hom_dim(pos[F(3, 4)], X), A.hom_dim(neg[F(3, 4)], X)
    names = ['S18', 'S25', 'S69', 'S74']
    span_pos, span_nine = [], []
    for mask in range(16):
        sel = [names[i] for i in range(4) if mask >> i & 1]
        b = []
        for s in sel:
            b = A.xor(b, O.switch7(s))
        X = A.deform(N7, b)
        assert A.is_mc(X)
        p = pr(X, +1)
        if p == (9, 31):
            span_pos.append('+'.join(sel) or '0')
        if p[0] == 9:
            span_nine.append('+'.join(sel) or '0')
    rep.check('smoothings with the pair (9, 31) for t > 0',
              ['S18 census', 'S18 opposite', 'S25 census', 'S25 opposite'], sorted(hits[+1]))
    rep.check('span elements with the pair (9, 31) for t > 0', ['S18', 'S25'], sorted(span_pos))
    rep.check('span elements reaching dimension nine for t > 0', ['S18', 'S25'], sorted(span_nine))
    rep.check('third closure for t > 0: opposite smoothings at S_18, S_25 against the census smoothings', (69, 69, 103, 103),
              (third[('S18', 'opposite')][0], third[('S25', 'opposite')][0], third[('S18', 'census')][0],
               third[('S25', 'census')][0]))
    rep.check('the same for t < 0 (Proposition 6.11)', (69, 69, 103, 103),
              (third[('S18', 'opposite')][1], third[('S25', 'opposite')][1], third[('S18', 'census')][1],
               third[('S25', 'census')][1]), 'expected')
    rep.check('smoothings with the pair (9, 31) for t < 0 (Propositions 6.10, 6.11)',
              ['S18 census', 'S18 opposite', 'S25 census', 'S25 opposite'], sorted(hits[-1]), 'expected')
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
