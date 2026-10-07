"""Remark 4.11 (the cases q >= 11): finite searches among reduced loop-type twisted complexes over F2.

N_q is the encoding of the piecewise-linear model (q = 5, 7, 11, 13) or, for q = 17, 19, the encoding of a numerical
trace of Smith's curve (data/traced_curves.txt; kwz_traced_curves.py checks that the two agree for q <= 13).  The
target is BN_q.  Each search is exhaustive within the class it states.
  1. Edit distance from N_q to BN_q: 2 for q = 5, 7; 4 for q = 11, 13; at least 6 for q = 17, 19 (the search
     excludes edit distances at most 5).
  2. Chains of algebraic switches (each Maurer--Cartan for the complex it deforms) followed by one replacement of an
     arrow: at q = 5, 7 one switch suffices; at q = 11 and 13 no chain of at most three switches ends one replacement
     away from BN_q; at q = 11 the first such chains have four switches, and they end at the sum of an arc and a
     closed curve other than alpha_11 + (c, J_2).
  3. Switch distance from N_q to alpha_q + (c, J_2): 1 at q = 5, 7; more than 5 at q = 11.
  4. The 35 004 complexes (6 331 up to strict isomorphism) obtained from N_11 by at most two algebraic switches
     (N_11, its 268 switches, and the 34 736 distinct complexes obtained by a second switch, among which N_11 recurs): none has alpha_11 as a summand, and each pairs with E_{-1/2} to a
     number other than dim I(P(-2,3,11)) = 13, or pairs with the earring of some slope of a panel of 25 slopes to more
     than the rank of reduced Khovanov homology of the closure.
Option --quick: parts 1 to 3 for q <= 13 without the four-switch search at q = 11 and the switch distance at q = 11;
part 4 only for the first 2000 objects.
"""
import argparse, sys, time
from collections import Counter
from fractions import Fraction as F

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O
import switches_common as S

PANEL = [F(1, 2), F(3, 4), F(-1, 2), F(-3, 4), F(1), F(-1), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(0), None,
         F(1, 4), F(-1, 4), F(3, 2), F(-3, 2), F(2), F(-2), F(1, 5), F(2, 5), F(3, 5), F(4, 5), F(5, 4), F(10, 29),
         F(-10, 29)]                          # paper slopes r of the panel


def N_of(q):
    if q in (5, 7):
        return O.N(q)[0]
    if q in (11, 13):
        return O.complexes('pl_models')[f'N{q}'][0]
    traced = O.complexes('traced_curves')
    return traced[sorted(k for k in traced if k.startswith(f'q{q}_'))[0]][0]


def chain_search(N0, T, mmax, log):
    """Smallest m <= mmax such that some chain of m algebraic switches from N0 ends one replacement away from T, and
    the objects Y at which the chains with that m end.  Meet in the middle: S_a(N0) against b switches from K(T)."""
    Nside, _ = S.layers(N0, 2)
    K = S.one_replacement_neighbours(T, S.label_multiset(N0))
    log(f'|S_1|, |S_2| = {sum(1 for v in Nside.values() if v == 1)}, {sum(1 for v in Nside.values() if v == 2)}; '
        f'|K(T)| = {len(K)}')
    for m in range(0, min(mmax, 2) + 1):
        hits = [Y for (Y, t, tp, corner, c) in K if Nside.get(c) == m]
        if hits:
            return m, hits
    for m in range(3, mmax + 1):
        hits = []
        for (Y, t, tp, corner, c) in K:
            found = False
            for b, Z in S.algebraic_switches(Y):
                cz = A.canon(Z)
                if m == 3 and Nside.get(cz) == 2:
                    found = True
                    break
                if m == 4:
                    for b2, W in S.algebraic_switches(Z):
                        if Nside.get(A.canon(W)) == 2:
                            found = True
                            break
                    if found:
                        break
            if found:
                hits.append(Y)
        if hits:
            return m, hits
    return None, []


def switch_distance_to(N0, Y0, bmax):
    """Least number of switches joining N0 and Y0 if it is at most 2 + bmax, else None."""
    SN, _ = S.layers(N0, 2)
    seen = {A.canon(Y0): 0}
    fr = [Y0]
    best = [SN[A.canon(Y0)]] if A.canon(Y0) in SN else []
    for d in range(1, bmax + 1):
        nxt = []
        for Z in fr:
            for b, W in S.algebraic_switches(Z):
                c = A.canon(W)
                if c in SN:
                    best.append(SN[c] + d)
                if c not in seen:
                    seen[c] = d
                    if d < bmax:
                        nxt.append(W)
        fr = nxt
    return min(best) if best else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args()
    t0 = time.time()
    rep = A.Report('Remark 4.11: the cases q >= 11')
    log = lambda s: rep.note(f'{s} ({time.time() - t0:.0f} s)')

    rep.section('1. edit distance from N_q to BN_q')
    exp = {5: 2, 7: 2, 11: 4, 13: 4}
    for q in (5, 7, 11, 13):
        r, lb, eds, _ = S.edit_distance(N_of(q), G.bn_smith(q), rmax=6, timecap=600)
        rep.check(f'q={q}: edit distance', exp[q], r)
    if not args.quick:
        for q in (17, 19):
            X = N_of(q)
            r, lb, eds, excluded = S.edit_distance(X, G.bn_smith(q), rmax=5, timecap=1800, all_edits=False)
            rep.check(f'q={q}: no edit with at most 5 replacements, i.e. edit distance at least 6', (None, 5),
                      (r, excluded))

    rep.section('2. chains of algebraic switches followed by one replacement')
    for q in (5, 7):
        m, hits = chain_search(N_of(q), G.bn_smith(q), 1, log)
        rep.check(f'q={q}: least number of switches', 1, m)
    for q in (11, 13):
        m, hits = chain_search(N_of(q), G.bn_smith(q), 3, log)
        rep.check(f'q={q}: no chain of at most three switches', None, m)
    if not args.quick:
        m, hits = chain_search(N_of(11), G.bn_smith(11), 4, log)
        rep.check('q=11: least number of switches', 4, m)
        shapes = Counter(tuple(sorted(('closed' if c else 'arc', len(s)) for s, l, c in A.walk(Y))) for Y in hits)
        rep.note(f'end objects of the four-switch chains (components): {dict(shapes)}')
        ref = A.union(G.alpha(11), G.c_curve(2))
        rep.check('q=11: each chain ends at the sum of an arc and one closed curve, not alpha_11 + (c, J_2)', True,
                  all(sum(1 for s, l, c in A.walk(Y) if c) == 1 and not A.strictly_isomorphic(Y, ref) for Y in hits))

    rep.section('3. switch distance from N_q to alpha_q + (c, J_2)')
    for q in (5, 7):
        rep.check(f'q={q}: switch distance', 1, switch_distance_to(N_of(q), A.union(G.alpha(q), G.c_curve(2)), 3))
    if not args.quick:
        rep.check('q=11: switch distance > 5 (none within 2 + 3)', None,
                  switch_distance_to(N_of(11), A.union(G.alpha(11), G.c_curve(2)), 3))

    rep.section('4. the objects obtained from N_11 by at most two algebraic switches')
    N11 = N_of(11)
    singles = S.algebraic_switches(N11)
    rep.check('q=11: algebraic switches of N_11', 268, len(singles))
    seen, two = set(), []
    for b, Y in singles:
        for b2, Y2 in S.algebraic_switches(Y):
            key = frozenset(Y2[1])
            if key in seen:
                continue
            seen.add(key)
            two.append(Y2)
    rep.check('q=11: distinct objects obtained by a second switch', 34736, len(two))
    objs = [N11] + [Y for b, Y in singles] + two
    rep.note(f'objects listed: {len(objs)} (N_11 recurs once among the second switches)')
    rep.check('q=11: distinct complexes', 35004, len({frozenset(X[1]) for X in objs}))
    rep.check('q=11: distinct up to strict isomorphism', 6331, len({A.canon(X) for X in objs}))
    al = G.alpha(11)
    rep.check('q=11: objects with alpha_11 as a summand', 0,
              sum(1 for X in objs for c in A.components(X) if len(c[0]) == len(al[0]) and A.strictly_isomorphic(c, al)))
    kh = O.khovanov_ranks()
    ranks = {r: kh[(O.slope_key(r), 11)] for r in PANEL}
    E = {r: G.earring(r) for r in PANEL}
    rep.check('q=11: BN_11 pairs to rk Kh~ at every panel slope', True,
              all(A.hom_dim(E[r], G.bn_smith(11)) == ranks[r] for r in PANEL), 'expected')
    if args.quick:
        objs = objs[:2000]
        rep.note('--quick: only the first 2000 objects are screened')
    dist = Counter()
    compatible = 0
    for t, X in enumerate(objs):
        v = A.hom_dim(E[F(-1, 2)], X)
        dist[v] += 1
        if v == ranks[F(-1, 2)] and all(A.hom_dim(E[r], X) <= ranks[r] for r in PANEL):
            compatible += 1
        if t and t % 5000 == 0:
            log(f'{t} objects screened')
    rep.note(f'pairings with E_-1/2: {dict(sorted(dist.items()))}')
    rep.check('q=11: rk Kh~(P(-2,3,11)) = dim I(P(-2,3,11))', 13, ranks[F(-1, 2)])
    rep.check('q=11: objects pairing to 13 with E_-1/2 and to at most rk Kh~ at every panel slope', 0, compatible)
    if not args.quick:
        rep.check('q=11: pairings with E_-1/2 (distribution)', {11: 246, 13: 7356, 15: 25464, 17: 1939},
                  dict(sorted(dist.items())), 'expected')
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
