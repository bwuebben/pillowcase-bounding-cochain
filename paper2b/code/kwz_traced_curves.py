"""Traced curves (Appendix A.4, "Traced curves"; Remark 4.11, first sentences; Appendix A.1).

The encodings in data/traced_curves.txt come from floating-point traces of Smith's perturbed curve (see the header of
that file); the statements below about them are corroboration, not proof.  Exact algebra is applied to the encoded
complexes.
  * For q <= 13 all traces at a given q have strictly isomorphic encodings, strictly isomorphic to N_q, the encoding
    of the piecewise-linear model L_{q,PL} (N_5 of (A.1), N_7 = N of Section 3.3, and the encodings of L_{11,PL},
    L_{13,PL} in data/pl_models.txt); at q = 11, 13 they have 47 and 55 generators.
  * The pairing with E_{-1/2} is q at q = 5, 7, q + 4 at q = 11, 13 and q + 8 at q = 17, 19, against
    rk Kh~(P(-2,3,q)) = q + 2; one switch raises it to q + 2 at q = 5, 7 (Theorem 4.1(a)).
"""
import sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def main():
    t0 = time.time()
    rep = A.Report('Traced curves (Appendix A.4, Remark 4.11)')
    traced = O.complexes('traced_curves')
    pl = O.complexes('pl_models')
    E = G.earring(F(-1, 2))
    kh = O.khovanov_ranks()
    model = {5: O.N5()[0], 7: O.N7()[0], 11: pl['N11'][0], 13: pl['N13'][0]}
    ngen = {11: 47, 13: 55}
    expected_pair = {5: 5, 7: 7, 11: 15, 13: 17, 17: 25, 19: 27}
    for q in (5, 7, 11, 13, 17, 19):
        names = sorted(k for k in traced if k.startswith(f'q{q}_'))
        objs = [traced[k][0] for k in names]
        rep.section(f'q = {q}: {len(names)} trace(s) ({", ".join(n[len(str(q)) + 2:] for n in names)})')
        rep.check(f'q={q}: the traces have strictly isomorphic encodings', True,
                  all(A.strictly_isomorphic(objs[0], X) for X in objs[1:]))
        if q in model:
            rep.check(f'q={q}: the encoding is strictly isomorphic to N_q', True, A.strictly_isomorphic(objs[0], model[q]))
        if q in ngen:
            rep.check(f'q={q}: number of generators', ngen[q], len(objs[0][0]))
        rep.check(f'q={q}: pairing with E_-1/2', expected_pair[q], A.hom_dim(E, objs[0]))
        rep.check(f'q={q}: rk Kh~(P(-2,3,{q}))', q + 2, kh[('-1/2', q)])
    for q in (5, 7):
        Nq, L = O.N(q)
        b = list(O.selected_b(q).values())[0]
        rep.check(f'q={q}: one switch raises the pairing to q + 2', q + 2, A.hom_dim(E, A.deform(Nq, b)))
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
