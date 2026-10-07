"""Finite searches among reduced loop-type twisted complexes (Remarks 4.2 and 4.11, Appendix A.4, "Searches").

* Algebraic switches (Remark 4.2): replace two arrows i -Y^k-> j and i' -Y^k-> j' of delta_X, with j and j' of the
  same idempotent (and i, i' of the same idempotent), by i -Y^k-> j' and i' -Y^k-> j, such that the result satisfies
  the Maurer--Cartan equation.
* Edit distance (Remark 4.11): the least number of arrows that must be deleted from delta_X, and the same number
  added, to reach a complex strictly isomorphic to a target arc T.  All arrows are non-identity elements, so every
  candidate is reduced and strict isomorphism is homotopy equivalence.  X + edit is strictly isomorphic to T with
  2r changed arrows if and only if the word of T is the concatenation of r + 1 pairwise disjoint intervals of the
  component words of X (arc intervals read in either direction, intervals of closed components cyclic); the deleted
  arrows are the arrows of X not internal to an interval and the added arrows join consecutive intervals.  A greedy
  longest-match factorization gives a lower bound (it is optimal for the factor-closed language of all intervals,
  ignoring disjointness), and a depth-first search over disjoint tilings decides the exact value.
* Chains of algebraic switches are enumerated breadth first up to strict isomorphism (canonical forms).
"""
import itertools
import time
from collections import Counter

import kwz_algebra as A


# ----------------------------------------------------------------------------------------------- switches

def algebraic_switches(X, require_mc=True):
    """All switches of X: list of (b, X + b), b = [old1, old2, new1, new2]."""
    g, a = X
    aset = set(a)
    out = []
    for t1, t2 in itertools.combinations(a, 2):
        (i, kd, k, j), (i2, kd2, k2, j2) = t1, t2
        if (kd, k) != (kd2, k2) or g[i] != g[i2] or g[j] != g[j2]:
            continue
        n1, n2 = (i, kd, k, j2), (i2, kd, k, j)
        if n1 in aset or n2 in aset or n1[0] == n1[3] or n2[0] == n2[3]:
            continue
        b = [t1, t2, n1, n2]
        Y = (list(g), A.xor(a, b))
        if require_mc and not A.is_mc(Y):
            continue
        out.append((b, Y))
    return out


def label_multiset(X):
    return Counter((X[0][t[0]], t[1], t[2]) for t in X[1])


def one_arc(X):
    return sum(1 for _, _, c in A.walk(X) if not c) == 1


def layers(X, depth, log=None):
    """{canonical form: least number of switches} for the objects reachable from X by chains of at most `depth`
    algebraic switches, every object of the chain Maurer--Cartan."""
    seen = {A.canon(X): 0}
    frontier = [X]
    for d in range(1, depth + 1):
        nxt = []
        for Z in frontier:
            for b, W in algebraic_switches(Z):
                c = A.canon(W)
                if c not in seen:
                    seen[c] = d
                    nxt.append(W)
        frontier = nxt
        if log:
            log(f'depth {d}: {len(nxt)} new objects')
    return seen, frontier


def one_replacement_neighbours(T, labels_needed):
    """K(T): the objects Y with one arc component, Maurer--Cartan, with the label multiset `labels_needed`, from which
    one replacement (delete one arrow, add one arrow) leads to T.  Returns a list of (Y, t, t', corner, canon(Y)),
    where Y = T - t + t' and corner records whether t is an R-chord at an end of the arc of Y."""
    g, a = T
    out = []
    seen = set()
    for t in a:
        rest = [x for x in a if x != t]
        deg = Counter()
        for (i, _, _, j) in rest:
            deg[i] += 1
            deg[j] += 1
        V1 = [v for v in range(len(g)) if deg[v] <= 1]
        for u, v in itertools.permutations(V1, 2):
            for (kd, k) in {(t[1], t[2])} | {(x[1], x[2]) for x in a}:
                if kd == 'D' and g[u] != g[v]:
                    continue
                if kd == 'S' and (g[u] == g[v]) != (k % 2 == 0):
                    continue
                tp = (u, kd, k, v)
                if tp == t:
                    continue
                Y = (g, rest + [tp])
                if label_multiset(Y) != labels_needed:
                    continue
                if not A.is_mc(Y) or not one_arc(Y):
                    continue
                c = A.canon(Y)
                if c in seen:
                    continue
                seen.add(c)
                ends = A.arc_ends(Y)
                corner = (t[1] == 'D' and g[t[0]] == 'b' and (t[0] in ends or t[3] in ends))
                out.append((Y, t, tp, corner, c))
    return out


# ----------------------------------------------------------------------------------------------- edit distance

def _sources(X):
    """Matching sources (generator sequence, labels, maximal interval length) for every component and orientation;
    closed components are doubled so that cyclic intervals become ordinary ones."""
    out = []
    for seq, labs, closed in A.walk(X):
        if not closed:
            out.append((seq, labs, len(seq)))
            out.append((seq[::-1], [A._flip(l) for l in reversed(labs)], len(seq)))
        else:
            n = len(seq)
            out.append((seq + seq, labs + labs, n))
            rs = seq[::-1]
            rl = [A._flip(labs[(n - 1 - i - 1) % n]) for i in range(n)]
            out.append((rs + rs, rl + rl, n))
    return out


class Tiler:
    """Tilings of the word of the arc T by disjoint intervals of the components of X."""

    def __init__(self, X, T):
        self.X = X
        self.g = X[0]
        (seqT, labsT, closed), = A.walk(T)
        assert not closed
        self.gT = [T[0][v] for v in seqT]
        self.aT = list(labsT)
        self.n = len(self.gT)
        assert len(self.g) == self.n
        self.src = _sources(X)
        self._pl = {}
        self._gb = {}
        self.nodes = 0

    def placements(self, p):
        if p in self._pl:
            return self._pl[p]
        out = []
        gT, aT, n, g = self.gT, self.aT, self.n, self.g
        for si, (seq, labs, cap) in enumerate(self.src):
            m = len(seq)
            for u in range(min(m, cap)):
                if g[seq[u]] != gT[p]:
                    continue
                L = 1
                while p + L < n and u + L < m and L < cap and labs[u + L - 1] == aT[p + L - 1] \
                        and g[seq[u + L]] == gT[p + L]:
                    L += 1
                out.append((L, si, u))
        out.sort(reverse=True)
        self._pl[p] = out
        return out

    def lower_bound(self, start=0):
        p, tiles = start, 0
        while p < self.n:
            pl = self.placements(p)
            if not pl:
                return 10 ** 6
            p += pl[0][0]
            tiles += 1
        return tiles - 1

    def _gbound(self, p):
        if p not in self._gb:
            self._gb[p] = self.lower_bound(p) + 1 if p < self.n else 0
        return self._gb[p]

    def tilings(self, r, cap=1000, deadline=None):
        """All disjoint tilings of T by exactly r + 1 intervals (up to cap).  Raises TimeoutError after deadline."""
        n = self.n
        sols = []
        failed = set()

        def dfs(p, used, left, acc):
            self.nodes += 1
            if deadline is not None and self.nodes % 100000 == 0 and time.time() > deadline:
                raise TimeoutError
            if len(sols) >= cap:
                return False
            if p == n:
                if left == 0:
                    sols.append(list(acc))
                    return True
                return False
            if left == 0 or self._gbound(p) > left:
                return False
            key = (p, used, left)
            if key in failed:
                return False
            found = False
            for Lmax, si, u in self.placements(p):
                seq = self.src[si][0]
                for L in range(Lmax, 0, -1):
                    if left == 1 and p + L != n:
                        continue
                    if p + L < n and self._gbound(p + L) > left - 1:
                        continue
                    m = 0
                    for v in seq[u:u + L]:
                        m |= 1 << v
                    if used & m:
                        continue
                    acc.append((p, L, si, u))
                    found |= dfs(p + L, used | m, left - 1, acc)
                    acc.pop()
            if not found:
                failed.add(key)
            return found

        dfs(0, 0, r + 1, [])
        return sols

    def to_edit(self, til):
        """A tiling -> (deleted arrows, added arrows), as arrows of X."""
        internal = set()
        firsts, lasts = [], []
        for (p, L, si, u) in til:
            seq, labs, cap = self.src[si]
            for i in range(L - 1):
                x, y, (kd, k, dr) = seq[u + i], seq[u + i + 1], labs[u + i]
                internal.add((x, kd, k, y) if dr > 0 else (y, kd, k, x))
            firsts.append(seq[u])
            lasts.append(seq[u + L - 1])
        deleted = [t for t in self.X[1] if t not in internal]
        added = []
        for t in range(len(til) - 1):
            p, L, si, u = til[t]
            kd, k, dr = self.aT[p + L - 1]
            x, y = lasts[t], firsts[t + 1]
            added.append((x, kd, k, y) if dr > 0 else (y, kd, k, x))
        return deleted, added


def edit_distance(X, T, rmax, timecap=None, all_edits=True, cap=1000):
    """Exact edit distance from X to the arc T, if it is at most rmax.

    Returns (r, lower bound, edits, excluded) where r is None if no edit with at most rmax replacements exists or the
    time cap was reached; `excluded` is the largest k such that every edit distance <= k has been excluded; edits are
    verified (Maurer--Cartan and strictly isomorphic to T)."""
    Tl = Tiler(X, T)
    lb = Tl.lower_bound()
    deadline = time.time() + timecap if timecap else None
    excluded = lb - 1
    if lb > rmax:
        return None, lb, [], excluded
    try:
        for r in range(lb, rmax + 1):
            sols = Tl.tilings(r, cap=cap if all_edits else 1, deadline=deadline)
            if sols:
                eds = {}
                for s in sols:
                    d, a = Tl.to_edit(s)
                    assert len(d) == len(a) == r
                    Y = (X[0], A.xor(X[1], d + a))
                    assert A.is_mc(Y) and A.strictly_isomorphic(Y, T), 'edit verification failed'
                    eds[(frozenset(d), frozenset(a))] = (d, a)
                return r, lb, list(eds.values()), r - 1
            excluded = r
    except TimeoutError:
        pass
    return None, lb, [], excluded
