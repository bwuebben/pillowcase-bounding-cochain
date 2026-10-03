#!/usr/bin/env python3
r"""Remark 4.3: a computer check of the base-case convention (requires SnapPy and spherogram).

The four-component link H_A u H_B u K u e of [HHK2, Fig. 18] is entered as a planar diagram code: K is the knot,
H_A and H_B are the two unknots whose surgeries produce the tangle complement, and e is the curve of the Conway
sphere enclosing the punctures a, b (the curve along which t_e twists), pushed into the trivial tangle.  Filling H_A
by -s/p and H_B by q/r with (p, q, r, s) = (3, 5, 2, -1) gives the exterior X of K^(0) u e.  The script checks:

  (1) X is isometric to the exterior Y of the closure of the positive braid (sigma_1 sigma_2)^5 together with an
      unknotted circle encircling two of its strands, by isometries that preserve orientation (the cusp maps have
      determinant +1) and carry meridians to meridians; SnapPy reports that they extend to the links.  An exact
      combinatorial isomorphism of canonical retriangulations is also sought;
  (2) as a further consistency check, for h = -8, ..., 1 the Dehn filling of the cusp e of X along the slope
      (1, h) gives a knot whose Alexander polynomial (by Fox calculus from SnapPy's presentation) is that of
      P(-2, 3, 5 - 2h) from (3.2); for h = 1 this is T(3,4) = P(-2,3,3).

Remark 4.3 records (1).  This check is not used in the proofs: the convention of Remark 4.3 is a premise of the paper.

With --verified (requires SnapPy inside SageMath, for interval arithmetic) the script gives a computer-assisted
proof that K^(0) is the positive torus knot T(3,5), relative to the diagram code PD_HHK below being a faithful
transcription of [HHK2, Fig. 18].  The triangulations of X and Y are produced from the diagrams by exact
combinatorial operations (Dehn filling included).  The checks are:

  (V0) surgery convention: the (-1, 1) filling of the circle cusp of Y is homeomorphic, preserving orientation and
       meridians, to the exterior of the closure of (sigma_1 sigma_2)^5 sigma_1^2.  By the Rolfsen twist, -1 surgery on
       that circle adds a right-handed full twist, so SnapPy's filling (a, b) is a/b surgery in the orientation in
       which the crossings of positive braids are right-handed (with the left-handed twist the closure would be T(3,4),
       a Seifert fibred knot);
  (V1) canonical_retriangulation(verified=True) verifies, by interval arithmetic, that the triangulations compared
       subdivide the Epstein-Penner canonical cell decompositions of X and Y.  Their combinatorial isomorphisms,
       computed exactly, are therefore all the isometries X -> Y.  The list is nonempty and every member has cusp maps
       of determinant +1 that fix the meridians.

An orientation-preserving homeomorphism of exteriors that takes meridians to meridians extends to an
orientation-preserving homeomorphism of S^3.  It carries K^(0) u e to the closure of (sigma_1 sigma_2)^5 together with
the circle, so K^(0) is the closure of a positive braid, the positive T(3,5), of signature -8.  Because the list in
(V1) is complete and contains no orientation-reversing map, X is not homeomorphic to Y by an orientation-reversing map:
the mirror of T(3,5) is excluded.
"""
import sys
import time
import warnings
from fractions import Fraction

warnings.filterwarnings('ignore')
VERIFIED = '--verified' in sys.argv

try:
    if VERIFIED:
        import sage.all  # noqa: F401  (SnapPy's verified routines need Sage)
    import snappy
    from spherogram.links.tangles import Tangle, Crossing, BraidTangle, IdentityBraid
except ImportError:
    print('SnapPy and spherogram are required (pip install snappy; with --verified, SnapPy inside Sage); skipped.')
    sys.exit(2)

import sympy as sp

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


# PD code of [HHK2, Fig. 18]: X[under_in, ., under_out, .] counterclockwise.
# Edges: K 1..12, H_B 13..18, H_A 19..22, e 23..26.
PD_HHK = [
    [12, 26, 1, 23],   # e (front) over K (strand from b to c)
    [25, 1, 26, 2],    # K over e (back)
    [5, 23, 6, 24],    # e (front) over K (strand from a to d)
    [24, 6, 25, 7],    # K over e (back)
    [2, 15, 3, 14],    # H_B over K
    [13, 4, 14, 3],    # K over H_B
    [11, 15, 12, 16],  # H_B over K
    [18, 8, 13, 9],    # K over H_B
    [7, 4, 8, 5],      # K over K
    [9, 20, 10, 19],   # H_A over K
    [20, 11, 21, 10],  # K over H_A
    [17, 22, 18, 19],  # H_A over H_B
    [21, 16, 22, 17],  # H_B over H_A
]


def hhk_two_cusp(p=3, q=5, r=2, s=-1):
    """Exterior of K u e after filling H_A by -s/p and H_B by q/r (cusp order of SnapPy: K, H_B, H_A, e)."""
    M = snappy.Link(PD_HHK).exterior()
    a, b = Fraction(-s, p), Fraction(q, r)
    M.dehn_fill([(0, 0), (b.numerator, b.denominator), (a.numerator, a.denominator), (0, 0)])
    return M.filled_triangulation()


def ring_tangle():
    """A (2,2)-tangle: two vertical strands and an unknotted circle around them."""
    A, B, C, D = Crossing('rA'), Crossing('rB'), Crossing('rC'), Crossing('rD')
    A[0] = C[3]; A[1] = B[3]; A[3] = C[0]
    C[2] = D[0]; B[1] = D[2]; B[0] = D[3]
    return Tangle(2, [A, B, C, D], [(C, 1), (D, 1), (A, 2), (B, 2)], 'Ring')


def braid_ring_link(word=(1, 2) * 5, n=3, ring_at=1):
    b = BraidTangle(list(word), n)
    R = IdentityBraid(ring_at - 1) | ring_tangle() | IdentityBraid(n - ring_at - 1)
    return (b * R).braid_closure()


t = sp.symbols('t')


def alexander(G):
    """Alexander polynomial of a knot group with a deficiency-one presentation, by Fox calculus."""
    gens = list(G.generators())
    rels = [[(ch.lower(), 1 if ch.islower() else -1) for ch in r] for r in G.relators()]
    E = sp.Matrix([[sum(e for g, e in r if g == x) for x in gens] for r in rels])
    v = E.nullspace()[0]
    v = v * sp.ilcm(*[sp.fraction(x)[1] for x in v])
    v = v / sp.igcd(*[int(x) for x in v])
    ex = {x: int(v[i]) for i, x in enumerate(gens)}
    J = []
    for r in rels:
        row = []
        for x in gens:
            acc, pref = 0, sp.Integer(1)
            for g, e in r:
                if g == x:
                    acc += pref if e > 0 else -pref * t ** (-ex[g])
                pref = pref * t ** (ex[g] * e)
            row.append(sp.expand(acc))
        J.append(row)
    J = sp.Matrix(J)
    n = len(gens)
    j = next(i for i in range(n) if ex[gens[i]] != 0)
    minor = J[:, [i for i in range(n) if i != j]].det() if n > 1 else sp.Integer(1)
    d = sp.factor(sp.cancel(minor * (t - 1) / (t ** abs(ex[gens[j]]) - 1)))
    num = sp.fraction(sp.together(d))[0]
    c = [int(x) for x in reversed(sp.Poly(sp.expand(num * t ** 400), t).all_coeffs())]
    while c and c[0] == 0:
        c.pop(0)
    while c and c[-1] == 0:
        c.pop()
    if c and c[-1] < 0:
        c = [-x for x in c]
    return c


def pretzel_alexander(q):
    """(3.2): 1 - t + sum_{j=3}^q (-1)^(j+1) t^j - t^(q+2) + t^(q+3); for q = 3 this is T(3,4), for q = 1 T(2,5)."""
    if q == 1:
        return [1, -1, 1, -1, 1]
    return [1, -1, 0] + [(-1) ** (j + 1) for j in range(3, q + 1)] + [0, -1, 1]


def exact_isomorphisms(A, B):
    """All combinatorial isomorphisms between the verified canonical retriangulations of A and B, with cusp data."""
    TA, TB = A.canonical_retriangulation(verified=True), B.canonical_retriangulation(verified=True)
    out = []
    for iso in TA.isomorphisms_to(TB):
        maps = [[list(map(int, row)) for row in m] for m in iso.cusp_maps()]
        dets = [m[0][0] * m[1][1] - m[0][1] * m[1][0] for m in maps]
        meridian = all(m[1][0] == 0 and abs(m[0][0]) == 1 for m in maps)
        out.append((iso.cusp_images(), maps, dets, meridian, iso.extends_to_link()))
    return out


def verified_main():
    t0 = time.time()
    L = snappy.Link(PD_HHK)
    check('component order of PD_HHK is K, H_B, H_A, e (12, 6, 4, 4 crossing visits)',
          [len(c) for c in L.link_components] == [12, 6, 4, 4])
    X, Y = hhk_two_cusp(), braid_ring_link().exterior()
    print(f'X: verified volume {X.volume(verified=True, bits_prec=100)};  Y: {Y.volume(verified=True, bits_prec=100)}')
    K = BraidTangle([1, 2] * 5, 3).braid_closure()
    check('the closure of (sigma_1 sigma_2)^5 is a positive diagram (10 crossings of sign +1) of signature -8',
          len(K.crossings) == 10 and all(c.sign == 1 for c in K.crossings) and K.signature() == -8)
    print('(V0) surgery convention: (-1,1) filling of the circle of Y versus closure of (sigma_1 sigma_2)^5 sigma_1^2')
    Z = Y.copy()
    Z.dehn_fill([(0, 0), (-1, 1)])
    W = BraidTangle([1, 2] * 5 + [1, 1], 3).braid_closure().exterior()
    res = exact_isomorphisms(Z.filled_triangulation(), W)
    for r in res:
        print(f'    cusp maps {r[1]}, determinants {r[2]}, meridian to meridian {r[3]}')
    check('-1 surgery adds a right-handed full twist (orientation- and meridian-preserving isomorphism exists)',
          len(res) > 0 and all(all(d == 1 for d in r[2]) and r[3] for r in res))
    print('(V1) all isometries X -> Y from verified canonical retriangulations')
    res = exact_isomorphisms(X, Y)
    for r in res:
        print(f'    cusp images {r[0]}, cusp maps {r[1]}, determinants {r[2]}, meridians to meridians {r[3]}, '
              f'extends to the links {r[4]}')
    check(f'the complete list ({len(res)} isometries) is nonempty; each preserves orientation, takes meridians to '
          'meridians and extends to the links', len(res) > 0 and all(all(d == 1 for d in r[2]) and r[3] and r[4]
                                                                     for r in res))
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


def main():
    if VERIFIED:
        return verified_main()
    t0 = time.time()
    X = hhk_two_cusp()
    Y = braid_ring_link().exterior()
    print(f'X: {X.num_cusps()} cusps, volume {float(X.volume()):.6f};  Y: {Y.num_cusps()} cusps, volume {float(Y.volume()):.6f}')
    print('(1) isometries X -> Y')
    K = BraidTangle([1, 2] * 5, 3).braid_closure()
    check('the braid closure in Y is positive: all 10 crossings of (sigma_1 sigma_2)^5 have sign +1',
          len(K.crossings) == 10 and all(c.sign == 1 for c in K.crossings))
    isos = X.is_isometric_to(Y, return_isometries=True)
    good = []
    for iso in isos:
        maps = [[list(map(int, row)) for row in m] for m in iso.cusp_maps()]
        dets = [m[0][0] * m[1][1] - m[0][1] * m[1][0] for m in maps]
        meridian = all(m[1][0] == 0 and abs(m[0][0]) == 1 for m in maps)
        print(f'    cusp images {iso.cusp_images()}, cusp maps {maps}, determinants {dets}, meridians to meridians '
              f'{meridian}, extends to the links {iso.extends_to_link()}')
        good.append(all(d == 1 for d in dets) and meridian and iso.extends_to_link())
    check('an isometry exists; every isometry found is orientation-preserving, carries meridians to meridians and '
          'extends to the links', len(isos) > 0 and all(good))
    try:
        found = None
        for i in range(20):
            A1, B1 = X.copy(), Y.copy()
            if i:
                A1.randomize()
                B1.randomize()
            TA, TB = A1._canonical_retriangulation(), B1._canonical_retriangulation()
            combos = TA.isomorphisms_to(TB)
            if combos:
                found = combos
                break
        if found is not None:
            print(f'    exact combinatorial isomorphisms of canonical retriangulations: {len(found)}')
            check('an exact combinatorial isomorphism of canonical retriangulations exists and extends to the links',
                  any(c.extends_to_link() for c in found))
        else:
            print('    no combinatorial isomorphism of canonical retriangulations found (not a failure of (1))')
    except Exception as exc:                               # private SnapPy interface; optional
        print(f'    canonical retriangulation not available in this SnapPy version ({type(exc).__name__})')
    print('(2) fillings of the cusp e along (1, h)')
    ok = True
    for h in range(-8, 2):
        M = X.copy()
        M.dehn_fill([(0, 0), (1, h)])
        a = alexander(M.filled_triangulation().fundamental_group())
        same = a == pretzel_alexander(5 - 2 * h)
        ok &= same
        print(f'    h = {h:+d}: Alexander polynomial {a} = that of P(-2,3,{5 - 2 * h}): {same}')
    check('the filling along (1, h) has the Alexander polynomial of P(-2,3,5-2h) for h = -8, ..., 1', ok)
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
