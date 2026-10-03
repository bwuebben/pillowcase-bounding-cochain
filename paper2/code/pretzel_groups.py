#!/usr/bin/env python3
r"""Lemmas 4.1, 4.2 and 5.7, and Appendix A.3(iii): the fundamental groups of the reglued knots.

For the tangle of Hedden, Herald and Kirk for T(3,n), with the parameters (r, s) of Section 5, the group of the tangle
complement is free on A, B and at eps = 0 the boundary meridians are

    a = A^(s+3) B^(n-r),  b = B^-r A^s,  c = B^-r a B^r,  d = a^-1 A^3 A^s B^-r A^-3 a          ([HHK2, (31)]),

so that (r, s) = (nu+1, -1) gives (4.2) for n = 3 nu + 2, and (r, s) = (2 nu + 1, -2) gives the words of Lemma 5.7(b)
for n = 3 nu + 1.  By Lemma 4.1 the group of the knot K^(h) obtained by regluing with t_e^h is <A, B | t_e^h(b) = c>,
t_e^h(b) = (ba)^h b (ba)^-h.  The script checks, in the free group (exact reduced words):

  (1) ba = cd for these words (Appendix A.3(iii)), and at h = 0 the relation b = c is A^3 B^5 = 1 for n = 5;
  (2) Lemma 4.2, recorded after its proof: for n in {5, 8, ..., 29} and 0 <= m <= 15, the relation t_e^-m(b) = c is
      equivalent to the cyclic word R_{n,m} = B^(nu+1) A C^m A^2 B^(2nu+1) C^-m, C = A B^nu, and under alpha = A,
      beta = B^-1 the relator of Clay and Watson (4.4) is a cyclic permutation of R_{n,m}; with t_e^+m in place of
      t_e^-m (m >= 1) the relator of (4.4) is never obtained;
  (3) Lemma 5.7(b), as a word identity: for n = 3 nu + 1 with 1 <= nu <= 15 and 1 <= m <= 15, alpha -> A B^(2nu+1),
      beta -> B carries the relator of Clay and Watson for T(3, 3nu+2; 2, m-1) to a cyclic permutation of the
      relator of <A, B | t_e^-m(b) = c> or of its inverse (the lemma is proved by hand; this is a check);
  (4) Appendix A.3(iii): for every n <= 26 prime to 3 (n >= 4, where the tangle is defined) and 0 <= m <= 6, the
      Alexander polynomial computed by Fox calculus from <A, B | t_e^-m(b) = c> equals that of the closure of
      (sigma_1 sigma_2)^n sigma_1^(2m), computed with the reduced Burau representation.

All computations are exact.
"""
import sys
import time

import sympy as sp

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


# ------------------------------------------------------------------------- free group words: lists of (letter, +-1)
def red(w):
    out = []
    for g in w:
        if out and out[-1][0] == g[0] and out[-1][1] == -g[1]:
            out.pop()
        else:
            out.append(g)
    return out


def inv(w):
    return [(g, -e) for g, e in reversed(w)]


def pw(w, n):
    return red(w * n) if n >= 0 else red(inv(w) * (-n))


def mul(*ws):
    out = []
    for w in ws:
        out = red(out + list(w))
    return out


def cyc(w):
    w = red(w)
    while len(w) >= 2 and w[0][0] == w[-1][0] and w[0][1] == -w[-1][1]:
        w = w[1:-1]
    return w


def cyc_eq(u, v, allow_inverse=True):
    """u and v agree up to cyclic permutation (and, if allowed, inversion)."""
    u, v = cyc(u), cyc(v)
    if len(u) != len(v):
        return False
    cands = [v] + ([cyc(inv(v))] if allow_inverse else [])
    return any(c[i:] + c[:i] == u for c in cands for i in range(len(c)))


A, B = [('A', 1)], [('B', 1)]


def rs(n):
    """The parameters (r, s) of the tangle for T(3, n) used in Sections 4 and 5."""
    return ((2 * n + 1) // 3, -2) if n % 3 == 1 else ((n + 1) // 3, -1)


def words(n):
    r, s = rs(n)
    a = mul(pw(A, s + 3), pw(B, n - r))
    b = mul(pw(B, -r), pw(A, s))
    c = mul(pw(B, -r), a, pw(B, r))
    d = mul(inv(a), pw(A, 3), pw(A, s), pw(B, -r), pw(A, -3), a)
    return a, b, c, d


def relator(n, h):
    """The relator t_e^h(b) c^-1 of Lemma 4.1."""
    a, b, c, d = words(n)
    ba = mul(b, a)
    return mul(pw(ba, h), b, pw(ba, -h), inv(c))


def cw_relator(nu, m, alpha, beta):
    """Clay and Watson's relator for T(3, 3nu+2; 2, m) [CW, Prop. 22], with alpha, beta replaced by words."""
    y = mul(pw(beta, -nu), alpha)
    lhs = mul(pw(alpha, 2), pw(y, m), alpha)
    rhs = mul(pw(beta, 2 * nu + 1), pw(y, m), pw(beta, nu + 1))
    return mul(lhs, inv(rhs))


# ------------------------------------------------------------------------- Alexander polynomials
t = sp.symbols('t')


def normalize(expr):
    n, d = sp.fraction(sp.together(expr))
    P = sp.Poly(sp.expand(n), t)
    assert sp.Poly(sp.expand(d), t).is_monomial
    c = P.all_coeffs()[::-1]
    while c and c[0] == 0:
        c = c[1:]
    while c and c[-1] == 0:
        c = c[:-1]
    if c and c[0] < 0:
        c = [-x for x in c]
    return tuple(int(x) for x in c)


def fox_alexander(rel):
    """Alexander polynomial of the knot group <A, B | rel> by Fox calculus."""
    nA = sum(e for g, e in rel if g == 'A')
    nB = sum(e for g, e in rel if g == 'B')
    from math import gcd
    g0 = gcd(nA, nB)
    ab = {'A': nB // g0, 'B': -nA // g0}               # abelianization A -> t^ab[A], B -> t^ab[B]

    def phi(w):
        return t ** sum(ab[x] * e for x, e in w)
    acc, pref = sp.Integer(0), []
    for x, e in rel:
        if x == 'A':
            acc += phi(pref) if e == 1 else -phi(pref + [(x, -1)])
        pref = pref + [(x, e)]
    return normalize(sp.cancel(sp.factor(sp.expand(acc)) * (t - 1) / (t ** ab['B'] - 1)))


def burau_alexander(n, m):
    s1 = sp.Matrix([[-t, 1], [0, 1]])
    s2 = sp.Matrix([[1, 0], [t, -t]])
    M = (s1 * s2) ** n * s1 ** (2 * m)
    return normalize(sp.cancel(sp.factor((sp.eye(2) - M).det()) / (1 + t + t ** 2)))


def main():
    t0 = time.time()
    print('(1) the boundary relation ba = cd (Appendix A.3(iii))')
    good = []
    for n in [n for n in range(4, 30) if n % 3]:
        a, b, c, d = words(n)
        good.append(mul(b, a) == mul(c, d))
    check('ba = cd for every n in [4, 29] prime to 3, with the parameters (r, s) above', all(good))
    check('n = 5, h = 0: b = c reads A^3 B^5 = 1', cyc_eq(relator(5, 0), mul(pw(A, 3), pw(B, 5))))

    print('(2) Lemma 4.2: n = 3 nu + 2 in {5, 8, ..., 29}, 0 <= m <= 15')
    for nu in range(1, 10):
        n = 3 * nu + 2
        Cw = mul(A, pw(B, nu))
        ok1 = ok2 = ok3 = True
        for m in range(0, 16):
            R = mul(pw(B, nu + 1), A, pw(Cw, m), pw(A, 2), pw(B, 2 * nu + 1), pw(Cw, -m))
            ok1 &= cyc_eq(relator(n, -m), R)
            ok2 &= cyc_eq(R, cw_relator(nu, m, A, inv(B)), allow_inverse=False)
            if m >= 1:
                ok3 &= not cyc_eq(relator(n, m), cw_relator(nu, m, A, inv(B)))
        check(f'n = {n:2d}: t_e^-m(b) = c  <=>  R_(n,m) = 1; (4.4) under alpha = A, beta = B^-1 is a cyclic '
              f'permutation of R_(n,m); t_e^+m never gives (4.4)', ok1 and ok2 and ok3)

    print('(3) Lemma 5.7(b) as a word identity: n = 3 nu + 1, 1 <= nu <= 15, 1 <= m <= 15')
    ok = True
    for nu in range(1, 16):
        n = 3 * nu + 1
        for m in range(1, 16):
            img = cw_relator(nu, m - 1, mul(A, pw(B, 2 * nu + 1)), B)
            ok &= cyc_eq(img, relator(n, -m))
    check('alpha -> A B^(2nu+1), beta -> B carries the relator of T(3,3nu+2;2,m-1) to that of K_(-m)(3,3nu+1)', ok)
    ctrl = any(cyc_eq(cw_relator(nu, m, mul(A, pw(B, 2 * nu + 1)), B), relator(3 * nu + 1, -m))
               for nu in range(1, 6) for m in range(1, 6))
    check('control: with m in place of m - 1 no match occurs (1 <= nu, m <= 5)', not ctrl)

    print('(4) Fox calculus against Burau: 4 <= n <= 26 prime to 3, 0 <= m <= 6 (Appendix A.3(iii))')
    bad, cnt = [], 0
    for n in [n for n in range(4, 27) if n % 3]:
        for m in range(0, 7):
            cnt += 1
            if fox_alexander(relator(n, -m)) != burau_alexander(n, m):
                bad.append((n, m))
    check(f'{cnt} pairs (n, m): Fox(<A,B | t_e^-m(b) = c>) = Burau Alexander polynomial of (s1 s2)^n s1^(2m); '
          f'mismatches: {bad if bad else "none"}', not bad)
    print('  example n = 5, m = 1 (P(-2,3,7)):', fox_alexander(relator(5, -1)))

    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
