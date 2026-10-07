"""The two implementations of the morphism homology (Lemma 3.5 and Appendix A.4, "Morphism homology").

The dimensions of Propositions 3.7 and 3.8 are computed with both methods of kwz_algebra: the mapping-cone
reduction of Lemma 3.5 (with its data e, tau, r and the free ranks f_L, f_M, f_R) and the F2[H]-module structure.
Also checked: the data (e, tau, r) printed after (3.3); the torsion statements of Remark 3.6 (for End(N) one has
f_M = f_R = 1, for the one-generator complex (iota_circ, 0) one has f_L = f_M = 1, and p_{X,i} = U_X for every pair
used); and the agreement of the two methods on every pairing of Section 4 that the other scripts use.
"""
import sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def main():
    t0 = time.time()
    rep = A.Report('Morphism homology: Lemma 3.5 and the F2[H] method (Propositions 3.7, 3.8; Appendix A.4)')
    N, NL = O.N7()
    E, _ = O.E()
    Es, _ = O.Estar()
    objs = {'N': N}
    for s in ('S18', 'S25', 'S69', 'S74'):
        objs[s] = A.deform(N, O.switch7(s))
    paper = {'N': (7, 25), 'S18': (9, 31), 'S25': (9, 31), 'S69': (9, 23), 'S74': (9, 25)}
    etr = {'N': (228, 251, 227), 'S74': (228, 251, 227), 'S18': (228, 251, 224), 'S25': (228, 251, 224),
           'S69': (228, 251, 228)}
    rep.section('Propositions 3.7 and 3.8: (dim with E, dim with E_*)')
    for name, X in objs.items():
        rep.check(f'{name}: Maurer--Cartan', True, A.is_mc(X))
        h = (A.hom_dim(E, X), A.hom_dim(Es, X))
        c1, c2 = A.hom_dim_cone(E, X), A.hom_dim_cone(Es, X)
        rep.check(f'{name}: dimensions by the F2[H] method', paper[name], h)
        rep.check(f'{name}: dimensions by Lemma 3.5', paper[name], (c1['dim'], c2['dim']))
        rep.check(f'{name}: (e, tau, r) for the pairing with E_*', etr[name], (c2['e'], c2['tau'], c2['r']))
        rep.check(f'{name}: torsion condition f_L = f_M = f_R = 0 and p = U (both pairings)', True,
                  all(c['f'] == {'L': 0, 'M': 0, 'R': 0} and c['m'] == 1 for c in (c1, c2)))
    rep.section('Remark 3.6')
    cN = A.hom_dim_cone(N, N)
    rep.check('End(N): (f_L, f_M, f_R)', (0, 1, 1), (cN['f']['L'], cN['f']['M'], cN['f']['R']))
    one = (['w'], [])
    c1 = A.hom_dim_cone(one, one)
    rep.check('End(iota_circ, 0): (f_L, f_M, f_R)', (1, 1, 0), (c1['f']['L'], c1['f']['M'], c1['f']['R']))
    rep.section('the two methods on the pairings of Section 4 (additional check)')
    pairs = []
    for q in (5, 7):
        Nq, L = O.N(q)
        for bname, b in O.selected_b(q).items():
            X = A.deform(Nq, b)
            Y = A.deform(Nq, b, O.kappa(q))
            for r in (F(-1, 2), F(-3, 4), F(-1), O.r_q(q), F(1, 3)):
                Er = G.earring(r)
                pairs += [(f'q={q} {bname} E_{r}', Er, X), (f'q={q} {bname}+kappa E_{r}', Er, Y)]
        pairs += [(f'q={q} (c,J2) vs BN_q', G.c_curve(2), G.bn_smith(q)), (f'q={q} c vs BN_q', G.c_curve(1), G.bn_smith(q))]
    agree = 0
    for name, X, Y in pairs:
        d1 = A.hom_dim(X, Y)
        d2 = A.hom_dim_cone(X, Y)['dim']
        agree += d1 == d2
        if d1 != d2:
            rep.note(f'disagreement at {name}: {d1} vs {d2}')
    rep.check('pairings on which the two methods agree', len(pairs), agree, 'expected')
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
