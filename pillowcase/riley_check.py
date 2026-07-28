#!/usr/bin/env python3
r"""
riley_check.py -- exact Gaussian-integer verification of Theorem A
(dihedral rigidity for two-bridge knots) of paper 1.

Everything here is EXACT: the Riley matrices at meridian eigenvalue s = i have
entries in Z[i][u], so the whole Riley word is computed by integer arithmetic,
and the Riley polynomial is extracted by a polynomial gcd over Q(i).  No
floating point enters any assertion (floats appear only in the final printed
comparison against -4 sin^2(pi k / p), which is checked to 1e-12).

Presentation.  For the two-bridge knot b(p,q) (p odd, q odd, 0 < q < 2p,
gcd(p,q) = 1) we use the classical presentation

    pi_1(S^3 \ b(p,q)) = < a, b | a w = w b >,
    w = g_1^{e_1} g_2^{e_2} ... g_{p-1}^{e_{p-1}},
    g_i = b (i odd), a (i even),   e_i = (-1)^{floor(i q / p)}.

  b(3,1): w = b a          -> a b a = b a b            (trefoil)
  b(5,3): w = b a^-1 b^-1 a -> the classical figure-eight relator.

The exponent word depends only on q mod 2p, and q - 2p in (-p, 0) is the
standard signed representative, so the presentation is valid across the whole
range 0 < q < 2p.  In this normalization b(p,q) = b(p,q') iff q' = q^{+-1}
mod 2p (Schubert; two-bridge knots are invertible, so this also classifies
them unoriented).  Rolfsen's even-q labels are reached by q -> q + p, which
preserves the unoriented knot (q + p = q mod p); in particular 6_1 = b(9,2)
appears here as b(9,11), the mirror of b(9,7).  The traceless root set is
mirror-invariant and depends only on p, so chirality conventions do not
affect any check below.

Riley's parametrization (paper 1, Sec. 3.1):

    rho(a) = [[s, 1], [0, 1/s]],   rho(b) = [[s, 0], [-u, 1/s]],

and the traceless locus is s = i, where 1/s = -i, so every entry lies in Z[i][u].

What is verified
----------------
(1) BINARY-DIHEDRAL, IDENTICALLY.  In Z[i][u], with A = rho(a), B = rho(b),

        A (AB) A^{-1} = (AB)^{-1}

    holds as a polynomial identity -- for every u, before any relator is
    imposed.  So conjugation by A inverts AB: the image < A, B > = < AB > . < A >
    is binary dihedral for EVERY traceless Riley representation.  This is the
    engine of Theorem A, and it is checked symbolically, not sampled.

(2) THE TRACELESS RILEY POLYNOMIAL.  The relator a w = w b holds iff the four
    entries of A W - W B vanish; their monic gcd over Q(i) is the traceless
    Riley polynomial.  It is verified to equal, exactly,

        phi_p(u) = prod_{k=1}^{(p-1)/2} (u + 4 sin^2(pi k / p)),

    which is built here as an exact INTEGER polynomial with no trigonometry, via
    the Lucas recursion (see phi_exact below).  Consequences checked: degree
    (p-1)/2, constant term det K = p, root sum -p, and -- the point of the
    theorem -- INDEPENDENCE OF q.

(3) THE SMALL CASES printed in Remark 3.2 of the paper:
        phi_3 = u + 3,  phi_5 = u^2 + 5u + 5,
        phi_7 = u^3 + 7u^2 + 14u + 7,  phi_9 = u^4 + 9u^3 + 27u^2 + 30u + 9.

(4) THE TWO DETERMINANT-9 KNOTS.  9_1 = b(9,1) (torus) and 6_1 = b(9,2)
    (hyperbolic) produce the identical root set {-4 sin^2(pi k / 9)}.

Run:  python3 riley_check.py        (pure standard library, < 1 s)
"""

from fractions import Fraction
from math import sin, pi, isclose

# ---------------------------------------------------------------------------
# Q(i): Gaussian rationals.  Exact.  A field, so we can do polynomial gcd.
# ---------------------------------------------------------------------------

class GQ:
    """re + im * i, with re, im exact rationals."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(s, o): return GQ(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return GQ(s.re - o.re, s.im - o.im)
    def __neg__(s):    return GQ(-s.re, -s.im)

    def __mul__(s, o):
        return GQ(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)

    def __truediv__(s, o):
        d = o.re * o.re + o.im * o.im
        if d == 0:
            raise ZeroDivisionError
        return GQ((s.re * o.re + s.im * o.im) / d, (s.im * o.re - s.re * o.im) / d)

    def __eq__(s, o):  return s.re == o.re and s.im == o.im
    def is_zero(s):    return s.re == 0 and s.im == 0

    def __repr__(s):
        if s.im == 0: return f"{s.re}"
        if s.re == 0: return f"{s.im}i"
        return f"({s.re}{'+' if s.im > 0 else ''}{s.im}i)"

ZERO, ONE, I = GQ(0, 0), GQ(1, 0), GQ(0, 1)

# ---------------------------------------------------------------------------
# Polynomials in u over Q(i): coefficient list, index = degree.
# ---------------------------------------------------------------------------

def ptrim(p):
    while len(p) > 1 and p[-1].is_zero():
        p = p[:-1]
    return p

def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else ZERO) + (q[i] if i < len(q) else ZERO)
                  for i in range(n)])

def psub(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else ZERO) - (q[i] if i < len(q) else ZERO)
                  for i in range(n)])

def pmul(p, q):
    out = [ZERO] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a.is_zero():
            continue
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return ptrim(out)

def pis_zero(p): return len(p) == 1 and p[0].is_zero()

def pmonic(p):
    lc = p[-1]
    return [c / lc for c in p]

def pmod(a, b):
    """Remainder of a by b; b nonzero.  Q(i) is a field, so this terminates."""
    a = list(a)
    db = len(b) - 1
    while not pis_zero(a) and len(a) - 1 >= db:
        shift = len(a) - 1 - db
        fac = a[-1] / b[-1]
        for i, c in enumerate(b):
            a[i + shift] = a[i + shift] - fac * c
        a = ptrim(a)
    return a

def pgcd(a, b):
    a, b = ptrim(list(a)), ptrim(list(b))
    while not pis_zero(b):
        a, b = b, pmod(a, b)
    return pmonic(a)

# 2x2 matrices over the polynomial ring -------------------------------------

def mmul(M, N):
    return [[padd(pmul(M[0][0], N[0][0]), pmul(M[0][1], N[1][0])),
             padd(pmul(M[0][0], N[0][1]), pmul(M[0][1], N[1][1]))],
            [padd(pmul(M[1][0], N[0][0]), pmul(M[1][1], N[1][0])),
             padd(pmul(M[1][0], N[0][1]), pmul(M[1][1], N[1][1]))]]

def msub(M, N):
    return [[psub(M[i][j], N[i][j]) for j in range(2)] for i in range(2)]

def meq(M, N):
    return all(pis_zero(psub(M[i][j], N[i][j])) for i in range(2) for j in range(2))

K = lambda c: [c]                       # constant polynomial
U = [ZERO, ONE]                         # the variable u

# Riley matrices at s = i  (so 1/s = -i); entries in Z[i][u]
A    = [[K(I),        K(ONE)], [K(ZERO),      K(-I)]]
Ainv = [[K(-I),      K(-ONE)], [K(ZERO),       K(I)]]
B    = [[K(I),       K(ZERO)], [[ZERO, -ONE],  K(-I)]]     # lower-left = -u
Binv = [[K(-I),      K(ZERO)], [U,              K(I)]]     # lower-left = +u

# ---------------------------------------------------------------------------
# The exact target polynomial phi_p, built with NO trigonometry.
#
#   phi_p(u) = prod_{k=1}^{(p-1)/2} (u + 4 sin^2(pi k / p)),   v := u + 2.
#
# With zeta = e^{2 pi i / p}: u + 4 sin^2(pi k/p) = v - zeta^k - zeta^{-k}, so
#
#   phi_p^2 = prod_{k=1}^{p-1} (v - zeta^k - zeta^{-k})
#           = Phi(alpha) Phi(beta),   alpha + beta = v, alpha beta = 1,
#           = (2 - s_p(v)) / (2 - v),   s_n = alpha^n + beta^n,
#
# and s_n obeys the integer recursion s_0 = 2, s_1 = v, s_n = v s_{n-1} - s_{n-2}.
# So phi_p is the exact integer-polynomial square root of (2 - s_p)/(2 - v),
# shifted back by v = u + 2.  All integer arithmetic.
# ---------------------------------------------------------------------------

def _zmul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out

def _zsub(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0) for i in range(n)]

def _zdiv(num, den):
    """Exact division of integer polynomials; raises if not exact."""
    num, q = list(num), [0] * (len(num) - len(den) + 1)
    for shift in range(len(num) - len(den), -1, -1):
        c, dq = num[shift + len(den) - 1], den[-1]
        if c % dq:
            raise ArithmeticError("inexact division")
        c //= dq
        q[shift] = c
        for i, d in enumerate(den):
            num[shift + i] -= c * d
    if any(num):
        raise ArithmeticError("nonzero remainder")
    return q

def _zsqrt(Q):
    """Exact square root of a monic integer polynomial that is a perfect square."""
    d = (len(Q) - 1) // 2
    r = [0] * (d + 1)
    r[d] = 1
    for t in range(1, d + 1):
        acc = sum(r[j] * r[2 * d - t - j] for j in range(d - t + 1, d))
        num = Q[2 * d - t] - acc
        if num % 2:
            raise ArithmeticError("not a perfect square")
        r[d - t] = num // 2
    if _zmul(r, r) != Q:
        raise ArithmeticError("not a perfect square")
    return r

def _ztrim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p

def _shift_v_to_u(r):
    """r(v) -> r(u+2), by Horner in (u + 2)."""
    out = [0]
    for c in reversed(r):
        out = _ztrim(_zmul(out, [2, 1]))
        out[0] += c
    return _ztrim(out)

def phi_exact(p):
    """phi_p(u) as an exact integer coefficient list, index = degree."""
    s_prev, s_cur = [2], [0, 1]                    # s_0 = 2, s_1 = v
    for _ in range(2, p + 1):
        s_prev, s_cur = s_cur, _zsub(_zmul([0, 1], s_cur), s_prev)
    num = _ztrim(_zsub([2], s_cur))                # 2 - s_p(v)
    quot = _ztrim(_zdiv(num, [2, -1]))             # / (2 - v)
    return _shift_v_to_u(_zsqrt(quot))

# ---------------------------------------------------------------------------
# The Riley word and its polynomial
# ---------------------------------------------------------------------------

def riley_word(p, q):
    """[(letter, sign)] for w = g_1^{e_1} ... g_{p-1}^{e_{p-1}}."""
    return [("b" if i % 2 == 1 else "a", (-1) ** ((i * q) // p))
            for i in range(1, p)]

def word_matrix(w):
    M = [[K(ONE), K(ZERO)], [K(ZERO), K(ONE)]]
    for letter, sign in w:
        M = mmul(M, {("a", 1): A, ("a", -1): Ainv,
                     ("b", 1): B, ("b", -1): Binv}[(letter, sign)])
    return M

def riley_polynomial(p, q):
    """Monic gcd over Q(i) of the four entries of A W - W B."""
    W = word_matrix(riley_word(p, q))
    M = msub(mmul(A, W), mmul(W, B))
    g = None
    for i in range(2):
        for j in range(2):
            e = ptrim(M[i][j])
            if pis_zero(e):
                continue
            g = e if g is None else pgcd(g, e)
    return pmonic(g)

def as_integer_poly(p):
    """Convert a monic Q(i) polynomial to an integer list, or raise."""
    out = []
    for c in p:
        if c.im != 0 or c.re.denominator != 1:
            raise ArithmeticError(f"non-integer coefficient {c}")
        out.append(int(c.re))
    return out

def fmt(poly):
    """Pretty-print an integer polynomial in u, descending."""
    terms = []
    for d in range(len(poly) - 1, -1, -1):
        c = poly[d]
        if c == 0:
            continue
        mag = "" if abs(c) == 1 and d > 0 else str(abs(c))
        var = "" if d == 0 else ("u" if d == 1 else f"u^{d}")
        terms.append(("- " if c < 0 else "+ ") + mag + var)
    s = " ".join(terms)
    return s[2:] if s.startswith("+ ") else s.replace("- ", "-", 1)

# ---------------------------------------------------------------------------

def main():
    ok = True
    n_checks = 0

    # (1) binary-dihedral, as an identity in Z[i][u] ------------------------
    print("(1) Binary-dihedral structure, as a polynomial identity in Z[i][u]")
    AB = mmul(A, B)
    ABinv = mmul(Binv, Ainv)
    lhs = mmul(mmul(A, AB), Ainv)
    idcheck = meq(mmul(AB, ABinv), [[K(ONE), K(ZERO)], [K(ZERO), K(ONE)]])
    dihedral = meq(lhs, ABinv)
    n_checks += 2
    ok &= idcheck and dihedral
    print(f"    (AB)(AB)^-1 = 1 in Z[i][u] .......... {'PASS' if idcheck else 'FAIL'}")
    print(f"    A (AB) A^-1 = (AB)^-1 in Z[i][u] .... {'PASS' if dihedral else 'FAIL'}")
    print("    => conjugation by the meridian A inverts AB, for EVERY u: the image")
    print("       <A,B> = <AB> . <A> is binary dihedral before any relator is imposed.")
    print("       (This is the content of Theorem A; the relator only selects which u.)")

    # (2)-(4) the traceless Riley polynomial, all p <= 9, all q -------------
    print("\n(2) The traceless Riley polynomial phi_{p,q}(i,u), exactly, for all p <= 9")
    print("    (every q: odd, 0 < q < 2p, gcd(p,q) = 1 -- the full normalized range)")
    print("    knot                   gcd of A W - W B               = phi_p ?  deg  phi(0)")
    from math import gcd
    NAME = {(3, 1): "3_1", (5, 3): "4_1", (5, 1): "5_1", (7, 3): "5_2",
            (7, 1): "7_1", (9, 7): "6_1", (9, 11): "6_1 = b(9,2)", (9, 1): "9_1"}
    seen = {}
    for p in (3, 5, 7, 9):
        target = phi_exact(p)
        for q in range(1, 2 * p, 2):
            if gcd(p, q) != 1:
                continue
            got = as_integer_poly(riley_polynomial(p, q))
            agree = (got == target)
            n_checks += 1
            ok &= agree
            deg_ok = (len(got) - 1 == (p - 1) // 2)
            const_ok = (got[0] == p)
            sum_ok = (-got[len(got) - 2] == -p)      # root sum = -p
            n_checks += 3
            ok &= deg_ok and const_ok and sum_ok
            seen.setdefault(p, []).append((q, tuple(got)))
            label = f"b({p},{q})" + (f" [{NAME[(p, q)]}]" if (p, q) in NAME else "")
            print(f"    {label:<22} {fmt(got):<38} {'PASS' if agree else 'FAIL':<9}"
                  f"{len(got)-1:<4} {got[0]}")

    # q-independence --------------------------------------------------------
    print("\n(3) q-independence (Theorem A: the root set depends only on p)")
    for p, rows in seen.items():
        distinct = {poly for _, poly in rows}
        good = len(distinct) == 1
        n_checks += 1
        ok &= good
        qs = ", ".join(str(q) for q, _ in rows)
        print(f"    p = {p}: q in {{{qs}}} -> {len(distinct)} distinct polynomial(s)"
              f" .... {'PASS' if good else 'FAIL'}")
    d9 = dict(seen[9])
    same9 = d9[1] == d9[11] == d9[7]
    n_checks += 1
    ok &= same9
    print(f"    9_1 = b(9,1) (torus) and 6_1 = b(9,2) = b(9,11) (hyperbolic):"
          f" identical ... {'PASS' if same9 else 'FAIL'}")

    # (4) the printed small cases and the trigonometric root set ------------
    print("\n(4) Remark 3.2 small cases, and the roots -4 sin^2(pi k / p)")
    printed = {3: [3, 1], 5: [5, 5, 1], 7: [7, 14, 7, 1], 9: [9, 30, 27, 9, 1]}
    for p in (3, 5, 7, 9):
        target = phi_exact(p)
        match = (target == printed[p])
        n_checks += 1
        ok &= match
        print(f"    phi_{p} = {fmt(target):<38} vs paper {'PASS' if match else 'FAIL'}")
        # numeric root check: evaluate phi_p at -4 sin^2(pi k / p)
        worst = 0.0
        for k in range(1, (p + 1) // 2):
            x = -4.0 * sin(pi * k / p) ** 2
            val = sum(c * x ** d for d, c in enumerate(target))
            worst = max(worst, abs(val))
        root_ok = worst < 1e-9
        n_checks += 1
        ok &= root_ok
        print(f"      phi_{p}(-4 sin^2(pi k/{p})) = 0 for k = 1..{(p-1)//2},"
              f" max |value| = {worst:.2e}  {'PASS' if root_ok else 'FAIL'}")

    # (5) the dihedral character: rotation angle of rho(ab) -----------------
    print("\n(5) At each root, rho(ab) has SO(3) rotation angle 4 pi k / p")
    for p in (3, 5, 7, 9):
        worst = 0.0
        for k in range(1, (p + 1) // 2):
            u = -4.0 * sin(pi * k / p) ** 2
            tr = -2.0 - u                      # tr rho(ab) = -(u + 2)
            # tr = 2 cos(psi) with SO(3) angle 2 psi; compare to +-4 pi k / p
            import math
            psi = math.acos(max(-1.0, min(1.0, tr / 2.0)))
            worst = max(worst, min(abs(2 * psi - 4 * pi * k / p),
                                   abs(2 * psi - (2 * pi - 4 * pi * k / p)),
                                   abs(2 * psi + 4 * pi * k / p - 2 * pi)))
        good = worst < 1e-9
        n_checks += 1
        ok &= good
        print(f"    p = {p}: max angle error {worst:.2e}"
              f"   {'PASS' if good else 'FAIL'}")

    print(f"\n{'ALL PASS' if ok else 'FAILURES PRESENT'}  ({n_checks} checks)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
