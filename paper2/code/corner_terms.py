"""Remark 7.2 (the corner term) and Appendix A.2 (corner terms and minimal edits).  Exact finite searches over F2.

1. One-arrow replacements.  For (N_q, b), b in {b_S18, b_S25} at q = 7 and b = b_5 at q = 5: every complex obtained
   by deleting one arrow and adding one arrow between two generators that meet at most one arrow, with every label
   that occurs in BN_q or in BN~(T_q) (both directions, all endpoint pairs).  The paper: 566 candidates for each b at
   q = 7 and 402 at q = 5; exactly two are strictly isomorphic to BN_q; in each the added arrow is labelled R and
   starts at the end of N_q over (0,0); kappa_7 is the one common to b_S18 and b_S25; the table of Appendix A.2.  The
   same search with BN~(T_q) in place of BN_q gives the terms for t > 0 listed in Appendix A.2.
2. Minimal edits.  The edit distance from N_q to BN_q is 2 for q = 5, 7, with two minimal edits in each case; at
   q = 7 both consist of arrows labelled R and touch the end 36.
3. Sums of switches (Remark 7.2).  N_7 has 110 algebraic switches and N_5 has 58; no sum of at most three of them
   (221 925 and 32 567 sums) is strictly isomorphic to BN_q.  Of the 45 distinct switches of census smoothings at
   q = 7, 41 are algebraic switches (65 of the 73 occurrences); the other four exchange the heads of an M-arrow and
   an M2-arrow, and no sum of at most three of the 45 is strictly isomorphic to BN_7.
"""
import itertools, os, sys, time
from collections import Counter

import kwz_algebra as A
import kwz_encode as G
import kwz_objects as O
import switches_common as S


def one_arrow_replacements(X, labels):
    """Yields (deleted arrow, added arrow, X') for every one-arrow replacement of X."""
    gens, arrows = X
    for dl in arrows:
        rest = [t for t in arrows if t != dl]
        deg = Counter()
        for (i, _, _, j) in rest:
            deg[i] += 1
            deg[j] += 1
        endp = [v for v in range(len(gens)) if deg[v] < 2]
        for u, v in itertools.permutations(endp, 2):
            for (kd, k) in labels:
                if kd == 'D' and gens[u] != gens[v]:
                    continue
                if kd == 'S' and (gens[u] == gens[v]) != (k % 2 == 0):
                    continue
                new = (u, kd, k, v)
                if new == dl:
                    continue
                yield dl, new, (gens, rest + [new])


def term(Nq, L, dl, new):
    """The element of End(N_q) as in Appendix A.2 (the first term is the deleted arrow)."""
    return f'{A.arrow_str(Nq, dl, L)} + {A.arrow_str(Nq, new, L)}'


def main():
    t0 = time.time()
    rep = A.Report('Remark 7.2 and Appendix A.2: corner terms, minimal edits, sums of switches')
    expected_ii = {(7, 'b_S18'): ['0 -R-> 19 + 36 -R-> 19', '22 -R-> 3 + 36 -R-> 3'],
                   (7, 'b_S25'): ['0 -R-> 19 + 36 -R-> 19', '4 -R-> 3 + 36 -R-> 3'],
                   (5, 'b_5'): ['12 -R-> 11 + 0 -R-> 11', '8 -R-> 15 + 0 -R-> 15']}
    expected_i = {7: ['2 -M-> 3 + 18 -M2-> 3', '20 -M-> 19 + 18 -M2-> 19'],
                  5: ['10 -M-> 11 + 22 -M2-> 11', '14 -M-> 15 + 22 -M2-> 15']}
    ncand = {7: 566, 5: 402}
    rep.section('1. one-arrow replacements of (N_q, b)')
    sols_ii = {}
    for q in (7, 5):
        Nq, L = O.N(q)
        BN, BNk = G.bn_smith(q), O.kht_bn(q)
        labels = sorted({(t[1], t[2]) for T in (BN, BNk) for t in T[1]})
        e0 = O.end_labels(q)['(0,0)']
        for bname, b in O.selected_b(q).items():
            X = A.deform(Nq, b)
            n = 0
            hit_ii, hit_i = [], []
            for dl, new, Y in one_arrow_replacements(X, labels):
                n += 1
                if not A.is_mc(Y) or not A.is_loop_type(Y):
                    continue
                if A.strictly_isomorphic(Y, BN):
                    hit_ii.append((dl, new))
                if A.strictly_isomorphic(Y, BNk):
                    hit_i.append((dl, new))
            rep.check(f'q={q}, {bname}: number of candidates', ncand[q], n)
            got = sorted(term(X, L, dl, new) for dl, new in hit_ii)
            rep.check(f'q={q}, {bname}: replacements giving BN_q', sorted(expected_ii[(q, bname)]), got)
            rep.check(f'q={q}, {bname}: in each, the added arrow is labelled R and starts at the end {e0} over (0,0)', True,
                      all(A.label_of(X, new) == 'R' and L[new[0]] == e0 for dl, new in hit_ii))
            sec = [dl for dl, new in hit_ii if dl not in Nq[1]]
            if q == 5 or bname == 'b_S18':
                rep.check(f'q={q}, {bname}: the second solution deletes an arrow of the smoothing',
                          ['8 -R-> 15'] if q == 5 else ['22 -R-> 3'], [A.arrow_str(X, t, L) for t in sec])
            sols_ii[(q, bname)] = got
            rep.check(f'q={q}, {bname}: replacements giving BN~(T_q) (t > 0)', sorted(expected_i[q]),
                      sorted(term(X, L, dl, new) for dl, new in hit_i))
    common = set(sols_ii[(7, 'b_S18')]) & set(sols_ii[(7, 'b_S25')])
    rep.check('q=7: the corner term common to b_S18 and b_S25 is kappa_7', ['0 -R-> 19 + 36 -R-> 19'], sorted(common))
    rep.check('q=5: kappa_5 is a solution', True, '12 -R-> 11 + 0 -R-> 11' in sols_ii[(5, 'b_5')])

    rep.section('2. minimal edits from N_q to BN_q')
    for q in (7, 5):
        Nq, L = O.N(q)
        r, lb, eds, _ = S.edit_distance(Nq, G.bn_smith(q), rmax=4)
        rep.check(f'q={q}: edit distance from N_q to BN_q', 2, r)
        rep.check(f'q={q}: number of minimal edits', 2, len(eds))
        for d, a in eds:
            rep.note('delete ' + ', '.join(A.arrow_str(Nq, t, L) for t in d) + '; add '
                     + ', '.join(A.arrow_str(Nq, t, L) for t in a))
        if q == 7:
            rep.check('q=7: both minimal edits consist of arrows labelled R and touch the end 36', True,
                      all(all(A.label_of(Nq, t) == 'R' for t in d + a) and any(36 in (L[t[0]], L[t[3]]) for t in a)
                          for d, a in eds))

    rep.section('3. sums of at most three algebraic switches')
    for q, nsw, nsum in ((7, 110, 221925), (5, 58, 32567)):
        Nq, L = O.N(q)
        BN = G.bn_smith(q)
        sw = S.algebraic_switches(Nq, require_mc=False)
        mc = sum(1 for b, Y in sw if A.is_mc(Y))
        rep.check(f'q={q}: number of algebraic switches of N_q', nsw, len(sw))
        rep.check(f'q={q}: each of them satisfies the Maurer--Cartan equation (number)', nsw, mc)
        cnt = hits = 0
        for r in (1, 2, 3):
            for combo in itertools.combinations(range(len(sw)), r):
                arr = Nq[1]
                for t in combo:
                    arr = A.xor(arr, sw[t][0])
                Y = (Nq[0], arr)
                cnt += 1
                if not A.is_loop_type(Y) or not A.is_mc(Y):
                    continue
                if len(Y[1]) == len(BN[1]) and A.strictly_isomorphic(Y, BN):
                    hits += 1
        rep.check(f'q={q}: number of sums of one, two or three switches', nsum, cnt)
        rep.check(f'q={q}: sums strictly isomorphic to BN_q', 0, hits)
    census = []
    for line in open(os.path.join(O.DATA, 'q7_census.txt'), encoding='utf-8'):
        if line.startswith('#') or not line.strip():
            continue
        name, kind, gp, sw = line.rstrip('\n').split('\t')
        if gp == 'yes':
            census.append((name, frozenset(O._terms(sw, O.N7_LABELS))))
    N7, L7 = O.N7()
    alg = {frozenset(b) for b, Y in S.algebraic_switches(N7)}
    rep.check('q=7: the generator-preserving census smoothings (number)', 73, len(census))
    distinct = {}
    for name, b in census:
        distinct.setdefault(b, []).append(name)
    other = sorted('/'.join(v) for b, v in distinct.items() if b not in alg)
    rep.check('q=7: census switches that are algebraic switches (occurrences, distinct)', (65, 41),
              (sum(b in alg for _, b in census), sum(b in alg for b in distinct)))
    rep.note('the other census switches exchange the heads of an M-arrow and an M2-arrow: ' + ', '.join(other))
    BN7 = G.bn_smith(7)
    keys = list(distinct)
    cnt = hits = 0
    for r in (1, 2, 3):
        for combo in itertools.combinations(keys, r):
            arr = N7[1]
            for b in combo:
                arr = A.xor(arr, list(b))
            Y = (N7[0], arr)
            cnt += 1
            if A.is_loop_type(Y) and A.is_mc(Y) and A.strictly_isomorphic(Y, BN7):
                hits += 1
    rep.check('q=7: sums of one, two or three of the 45 distinct census switches', 15225, cnt, 'expected')
    rep.check('q=7: such sums strictly isomorphic to BN_7', 0, hits, 'expected')
    return rep.finish(t0)


if __name__ == '__main__':
    sys.exit(main())
