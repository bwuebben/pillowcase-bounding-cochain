"""Example 3.19: two deformations of N = N_7 with the same two closure dimensions that a third curve separates.
Exact computation over F2.

  * (N, b_S18) and (N, b_S18 + b_S35) are Maurer--Cartan: S35 is a connector--inherited orbit, and b_S18 + b_S35
    is one of the 41 two-switch sums of Proposition 3.12 with the pair (9, 31);
  * both have pairing dimensions 9 with E = E_{-1/2} and 31 with E_* = E_{-3/4};
  * the figure-eight Z around the line of slope 1/3 through (0,0), which also passes through (pi,pi), pairs to 9
    with (N, b_S18) and to 11 with (N, b_S18 + b_S35), so the two objects are not homotopy equivalent;
  * the figure-eight around the parallel line through (pi,0) and (0,pi), the line of Q_{1/3}, pairs to 9 with both.
"""
import os
import sys

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))


def census():
    rows = {}
    with open(os.path.join(HERE, 'data', 'q7_census.txt')) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            name, kind, gp, arrows = line.rstrip('\n').split('\t')
            rows[name] = (kind, gp, arrows)
    return rows


def main():
    rep = A.Report('Example 3.19: a third test curve separates two objects with the same closure dimensions')
    N, L7 = O.N7()
    idx = {v: t for t, v in enumerate(L7)}
    E, _ = O.E()
    Es, _ = O.Estar()
    cen = census()
    rep.check('S35 is a connector--inherited orbit (connector--main in the data file) with a generator-preserving '
              'smoothing', ('connector--main', 'yes'), cen['S35'][:2])
    b18 = A.parse_arrows(cen['S18'][2], idx)
    b35 = A.parse_arrows(cen['S35'][2], idx)
    X1 = A.deform(N, b18)
    X2 = A.deform(N, b18, b35)
    rep.check('(N, b_S18 + b_S35) satisfies the Maurer--Cartan equation', True, A.is_mc(X2))
    rep.check('pairing dimensions of (N, b_S18) with E and E_*', (9, 31), (A.hom_dim(E, X1), A.hom_dim(Es, X1)))
    rep.check('pairing dimensions of (N, b_S18 + b_S35) with E and E_*', (9, 31), (A.hom_dim(E, X2), A.hom_dim(Es, X2)))
    w = G.direction(F(1, 3))
    Z = G.encode(G.figure8_points((0, 0), w), closed=True)
    rep.check('the line of Z joins the corners (pi,pi) and (0,0)', ('(pi,pi)', '(0,0)'),
              (G.corner_of((0, 0)), G.corner_of(w)))
    rep.check('pairing dimensions of Z with (N, b_S18) and (N, b_S18 + b_S35)', (9, 11),
              (A.hom_dim(Z, X1), A.hom_dim(Z, X2)))
    Zp = G.encode(G.figure8_points((1, 0), w), closed=True)
    rep.check('the parallel line joins (pi,0) and (0,pi)', ('(pi,0)', '(0,pi)'),
              (G.corner_of((1, 0)), G.corner_of((1 + w[0], w[1]))), source='expected')
    rep.check('the figure-eight around it pairs equally with both objects', (9, 9),
              (A.hom_dim(Zp, X1), A.hom_dim(Zp, X2)), source='expected')
    return rep.finish()


if __name__ == '__main__':
    sys.exit(main())
