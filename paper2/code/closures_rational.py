"""Rational closures: equation (7.1), Corollary 7.6 and the paragraph after it.  Exact computations over F2, with the
ranks of reduced Khovanov homology read from data/khovanov_ranks.txt (computed with Khoca 1.5).

For q = 5 (340 slopes) and q = 7 (342 slopes) of the slope set of Corollary 7.6:
  * (7.1): dim HF(E_r, BN_q) = rk Kh~(K_r(q); F2) at every slope (equivalently, by Phi, KWZ's pairing theorem
    dim HF(Kh~(Q_{-r}), BN~(T_q)) = rk Kh~(K_r(q)) in their chart, which is also checked);
  * Corollary 7.6: dim HF(E_r, (N_q, b)) = dim HF(E_r, (N_q, b + kappa_q)) = rk Kh~ at every slope, with the table
    at E_{-1} (parallel to c) and E_{r_q} (parallel to alpha_q), and the undeformed values 7, 28 (q = 5) and 11, 44
    (q = 7);
  * the slope set has 108 slopes with odd numerator and odd denominator, for which the arc of E_r joins (pi,0) to
    (0,pi); N_q pairs to less than rk Kh~ at 40 (q = 5) and 43 (q = 7) slopes and never to more; the common value is
    q + 2 at r = -1/2 and 31 at r = -3/4, q = 7.

Options:
  --khoca      recompute every rank with Khoca and SnapPy's spherogram (pip packages khoca, snappy) and compare
               with the data file;
  --sample N   with --khoca, recompute only N randomly chosen ranks per q.
"""
import argparse, random, sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def khoca_rank(r, q, _calc=[]):
    """rk Kh~(num(Q_r + Q_{1/3} + Q_{1/q}); F2) with Khoca (r = None: infinity)."""
    import warnings
    warnings.filterwarnings('ignore')
    from khoca import InteractiveCalculator
    from spherogram import RationalTangle
    if not _calc:
        _calc.append(InteractiveCalculator(2))
    Rx = RationalTangle(1, 0) if r is None else RationalTangle(r.numerator, r.denominator)
    L = (Rx + RationalTangle(1, 3) + RationalTangle(1, q)).numerator_closure()
    res = _calc[0]([list(c) for c in L.PD_code()])
    return sum(t[3] for t in res[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--khoca', action='store_true')
    ap.add_argument('--sample', type=int, default=0)
    args = ap.parse_args()
    t0 = time.time()
    rep = A.Report('Rational closures: (7.1), Corollary 7.6')
    kh = O.khovanov_ranks()
    nslopes = {5: 340, 7: 342}
    below = {5: 40, 7: 43}
    table = {5: {'-1': (7, 7, 7), 'rq': (30, 30, 30), 'N': (7, 28)},
             7: {'-1': (11, 11, 11), 'rq': (46, 46, 46), 'N': (11, 44)}}
    for q in (5, 7):
        rep.section(f'q = {q}')
        S = O.slope_set(q)
        rep.check('number of slopes', nslopes[q], len(S))
        odd = [r for r in S if G.parity(r) == 'odd/odd']
        rep.check('slopes with odd numerator and odd denominator', 108, len(odd))
        rep.check('for these the arc of E_r joins (pi,0) and (0,pi)', True,
                  all({G.corner_of(P) for P in (G.line_base(G.direction(-r), 'P1'),)} |
                      {G.corner_of((G.line_base(G.direction(-r), 'P1')[0] + G.direction(-r)[0],
                                    G.line_base(G.direction(-r), 'P1')[1] + G.direction(-r)[1]))}
                      == {'(pi,0)', '(0,pi)'} for r in odd))
        Nq, L = O.N(q)
        b = list(O.selected_b(q).values())
        Xs = [A.deform(Nq, bb) for bb in b]
        Ys = [A.deform(Nq, bb, O.kappa(q)) for bb in b]
        BN, BNk = G.bn_smith(q), O.kht_bn(q)
        rows = {}
        for r in S:
            key = O.slope_key(r)
            Er = G.earring(r)
            rows[key] = {'kh': kh[(key, q)], 'BN': A.hom_dim(Er, BN),
                         'BNkwz': A.hom_dim(G.kh_rational(None if r is None else -r), BNk),
                         'X': [A.hom_dim(Er, X) for X in Xs], 'Y': [A.hom_dim(Er, Y) for Y in Ys],
                         'N': A.hom_dim(Er, Nq)}
        rep.check('(7.1): dim HF(E_r, BN_q) = rk Kh~ (number of slopes)', nslopes[q],
                  sum(v['BN'] == v['kh'] for v in rows.values()))
        rep.check("KWZ's pairing theorem in their chart: dim HF(Kh~(Q_-r), BN~(T_q)) = rk Kh~ (number of slopes)",
                  nslopes[q], sum(v['BNkwz'] == v['kh'] for v in rows.values()), 'expected')
        rep.check('Corollary 7.6: the three numbers agree (number of slopes)', nslopes[q],
                  sum(all(x == v['kh'] for x in v['X'] + v['Y']) for v in rows.values()))
        e1, eq_ = rows['-1'], rows[str(O.r_q(q))]
        rep.check('E_{-1}: (N_q, b), (N_q, b + kappa_q), rk Kh~', table[q]['-1'], (e1['X'][0], e1['Y'][0], e1['kh']))
        rep.check(f'E_(r_q), r_q = {O.r_q(q)}: (N_q, b), (N_q, b + kappa_q), rk Kh~', table[q]['rq'],
                  (eq_['X'][0], eq_['Y'][0], eq_['kh']))
        rep.check('the undeformed N_q at E_{-1} and E_{r_q}', table[q]['N'], (e1['N'], eq_['N']))
        par = lambda u, v: u[0] * v[1] - u[1] * v[0] == 0
        rep.check('E_{-1} is parallel to c and E_{r_q} to alpha_q', True,
                  par(G.direction(F(1)), G.CDIR) and par(G.direction(-O.r_q(q)), (1 - q, 3 * q - 4)))
        rep.check('N_q pairs to less than rk Kh~ (number of slopes)', below[q],
                  sum(v['N'] < v['kh'] for v in rows.values()))
        rep.check('N_q pairs to more than rk Kh~ (number of slopes)', 0, sum(v['N'] > v['kh'] for v in rows.values()))
        rep.check('the common value at r = -1/2', q + 2, rows['-1/2']['kh'])
        if q == 7:
            rep.check('the common value at r = -3/4, q = 7', 31, rows['-3/4']['kh'])
        if args.khoca:
            keys = [r for r in S]
            if args.sample:
                keys = random.Random(q).sample(keys, min(args.sample, len(keys)))
            bad = [(O.slope_key(r), kh[(O.slope_key(r), q)], khoca_rank(r, q)) for r in keys]
            bad = [x for x in bad if x[1] != x[2]]
            rep.check(f'Khoca recomputation of {len(keys)} ranks: disagreements', [], bad, 'expected')
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
