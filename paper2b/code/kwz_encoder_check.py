"""The encoder of curves as twisted complexes, against Zibrowius's program kht++ (Appendix A.4 and Section 4.1).

1. For each of the 216 rational tangles Q_s in data/kht_rational.txt, BN~(Q_s) computed by kht++ is strictly
   isomorphic to the encoding of the straight lattice arc of KWZ slope s, and Kh~(Q_s) to the encoding of the
   figure-eight around it (Appendix A.4, "Encodings").
2. In Smith's coordinates, Phi(Kh~(Q_{-r})) is the earring curve E_r of Q^_r (Section 4.1): the encodings agree at
   49 slopes of the slope set of Corollary 4.6 (every seventh slope); we also check all 342 slopes.
3. The encodings of E_{-1/2} and E_{-3/4} are the complexes E and E_* printed in (3.2) and (3.3).
"""
import sys, time
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O


def kwz_slope(s):
    return None if s in ('∞', 'inf') else F(s)


def main():
    t0 = time.time()
    rep = A.Report('Encoder against kht++ (Appendix A.4, Section 4.1)')
    rep.section('1. rational tangles computed by kht++')
    data = O.read_kht_rational()
    ok_bn = ok_kh = 0
    bad = []
    for s, (bn, kh) in data.items():
        x = kwz_slope(s)
        b = A.strictly_isomorphic(bn, G.bn_rational(x))
        k = A.strictly_isomorphic(kh, G.kh_rational(x))
        ok_bn += b
        ok_kh += k
        if not (b and k):
            bad.append(s)
    rep.check('number of rational tangles compared', 216, len(data))
    rep.check('BN~(Q_s) is the straight lattice arc (number of slopes)', 216, ok_bn)
    rep.check('Kh~(Q_s) is the figure-eight around it (number of slopes)', 216, ok_kh)
    if bad:
        rep.note(f'slopes that fail: {bad}')

    rep.section('2. Phi(Kh~(Q_{-r})) and E_r in Smith\'s coordinates')
    S = O.slope_set(7)
    sample = S[::7]
    agree = sum(A.strictly_isomorphic(G.encode(G.shift(G.kh_rational_points(None if r is None else -r)), closed=True),
                                      G.earring(r)) for r in sample)
    rep.check('slopes in the sample', 49, len(sample))
    rep.check('Phi(Kh~(Q_{-r})) strictly isomorphic to E_r (number of sampled slopes)', 49, agree)
    allagree = sum(A.strictly_isomorphic(G.encode(G.shift(G.kh_rational_points(None if r is None else -r)), closed=True),
                                         G.earring(r)) for r in S)
    rep.check('the same at every slope of the set for q = 7', len(S), allagree, 'expected')

    rep.section('3. the printed earring complexes')
    E, _ = O.E()
    Es, _ = O.Estar()
    rep.check('E_{-1/2} strictly isomorphic to E of (3.2)', True, A.strictly_isomorphic(G.earring(F(-1, 2)), E))
    rep.check('E_{-3/4} strictly isomorphic to E_* of (3.3)', True, A.strictly_isomorphic(G.earring(F(-3, 4)), Es))
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
