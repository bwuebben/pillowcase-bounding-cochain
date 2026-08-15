#!/usr/bin/env python3
r"""
torus_characters.py -- verification of the (3,n)-torus knot traceless character
count (Proposition "Count" of paper 1, Sec. 4.2) and of the dihedral dichotomy
(Theorem B), by direct enumeration of the representation arcs.

Setup (paper 1, Sec. 4.1).  pi_1(S^3 \ T(3,n)) = < x, y | x^3 = y^n >, gcd(3,n)=1.
The element z = x^3 = y^n is central, so an irreducible SU(2) representation has

    rho(x) = exp(alpha u^),  alpha = pi l1 / 3,  l1 in {1,2},
    rho(y) = exp(beta  v^),  beta  = pi l2 / n,  0 < l2 < n,

with the central-sign match l1 = l2 (mod 2), and the angle theta between the
axes u^, v^ as the only remaining modulus: each admissible (l1,l2) is a
one-parameter arc of irreducible characters (Klassen).

The meridian is mu = x^a y^b with a n + 3 b = 1, and

    (1/2) tr rho(mu) = cos(a alpha) cos(b beta) - sin(a alpha) sin(b beta) cos(theta),

affine in cos(theta), so each arc carries exactly ONE traceless character, at

    cos(theta) = cot(a alpha) cot(b beta),                                  (*)

provided the right-hand side lies in (-1, 1); otherwise the arc carries none.

What is verified
----------------
(1) THE TRACE FORMULA.  For a dense sample of (l1, l2, theta), the closed form
    above is checked against a direct quaternion evaluation of rho(x)^a rho(y)^b.
    This is what makes (*) usable; it is checked, not assumed.

(2) THE COUNT N(3,n), by direct enumeration of the admissible arcs, against the
    closed formula of the Proposition, for every n <= 25 with gcd(3,n) = 1, and
    against the printed values N = 1,3,4,4,5,7,8,8 for n = 2,4,5,7,8,10,11,13.

(3) EACH TRACELESS REPRESENTATION IS GENUINE: for every admissible arc the
    representation is rebuilt as unit quaternions and checked to satisfy
    rho(x)^3 = rho(y)^n (central match), rho(z) = (-1)^{l1}, tr rho(mu) = 0,
    and irreducibility (theta strictly interior, so the axes are independent).

    The central sign is (-1)^{l1}, as the paper states (Sec. 4.2): the l1 = 1
    sheet has rho(x)^3 = rho(y)^n = -1, the l1 = 2 sheet has +1, and for n odd
    the traceless characters split evenly between the two sheets, exchanged by
    the symmetry (l1,l2) -> (3-l1, n-l2), which preserves (*).  The two sheets
    are reported separately below.

(4) THEOREM B (dihedral dichotomy).  rho(x) has trace 2 cos(pi l1 / 3) = +-1,
    never 0, so x never goes to a reflection; rho(y) is traceless iff l2 = n/2,
    possible only for n even.  Hence the number of irreducible traceless
    DIHEDRAL characters is 1 for n even and 0 for n odd, which is checked here
    to equal (det T(3,n) - 1)/2 with det = 3 (n even), 1 (n odd).

(5) N(3,n) = 2a FOR ALL ODD n <= 43, where a = -sigma(T(3,n))/4 is the number
    of irreducible flat SU(2) connections on the double branched cover
    Sigma(2,3,n), counted by the spherical triangle inequalities.  This is the
    relation used in Sec. 5 to pass from the character count to the signature,
    and n <= 43 matches the range claimed in the paper.  N is recomputed here
    by direct arc enumeration; the triangle count is imported from
    fs_gradings.py, so the two scripts are checked against each other rather
    than each against itself.

Run:  python3 torus_characters.py     (pure standard library, < 1 s)
"""

import math
from math import pi, cos, sin, gcd

TOL = 1e-11

# --- unit quaternions -------------------------------------------------------

def qmul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2)

def qinv(q):
    w, x, y, z = q
    return (w, -x, -y, -z)

def qpow(q, k):
    r = (1.0, 0.0, 0.0, 0.0)
    base = q if k >= 0 else qinv(q)
    for _ in range(abs(k)):
        r = qmul(r, base)
    return r

def qexp(angle, axis):
    """exp(angle * axis) for a unit imaginary axis."""
    ax, ay, az = axis
    s = sin(angle)
    return (cos(angle), s * ax, s * ay, s * az)

# --- the meridian -----------------------------------------------------------

def meridian_exponents(n):
    """(a, b) with a n + 3 b = 1."""
    for a in range(-3, 4):
        if (1 - a * n) % 3 == 0:
            return a, (1 - a * n) // 3
    raise ValueError(f"no meridian exponents for n = {n}")

def arcs(n):
    """All admissible (l1, l2): l1 in {1,2}, 0 < l2 < n, l1 = l2 mod 2."""
    return [(l1, l2) for l1 in (1, 2) for l2 in range(1, n)
            if (l1 - l2) % 2 == 0]

def cos_theta(n, l1, l2):
    """The right-hand side of (*), or None where a cotangent blows up."""
    a, b = meridian_exponents(n)
    alpha, beta = pi * l1 / 3.0, pi * l2 / n
    sa, sb = sin(a * alpha), sin(b * beta)
    if abs(sa) < 1e-14 or abs(sb) < 1e-14:
        return None
    return (cos(a * alpha) / sa) * (cos(b * beta) / sb)

def traceless_characters(n):
    """[(l1, l2, cos_theta)] for the arcs that carry a traceless character."""
    out = []
    for l1, l2 in arcs(n):
        c = cos_theta(n, l1, l2)
        if c is not None and -1.0 + 1e-12 < c < 1.0 - 1e-12:
            out.append((l1, l2, c))
    return out

def N_formula(n):
    """The closed form of the Proposition."""
    r = n % 6
    if r == 1: return 2 * (n - 1) // 3
    if r == 2: return (2 * n - 1) // 3
    if r == 4: return (2 * n + 1) // 3
    if r == 5: return 2 * (n + 1) // 3
    raise ValueError(f"gcd(3,{n}) != 1")

# ---------------------------------------------------------------------------

def main():
    ok = True
    n_checks = 0

    # (1) the closed trace formula vs direct quaternion evaluation ----------
    print("(1) The trace formula  (1/2)tr rho(mu) = cos(a al)cos(b be)"
          " - sin(a al)sin(b be)cos(th)")
    worst = 0.0
    for n in (5, 7, 8, 11, 13):
        a, b = meridian_exponents(n)
        for l1, l2 in arcs(n):
            alpha, beta = pi * l1 / 3.0, pi * l2 / n
            for j in range(1, 12):                       # sample theta in (0,pi)
                th = pi * j / 12.0
                X = qexp(alpha, (0.0, 0.0, 1.0))
                Y = qexp(beta, (sin(th), 0.0, cos(th)))
                MU = qmul(qpow(X, a), qpow(Y, b))
                closed = (cos(a * alpha) * cos(b * beta)
                          - sin(a * alpha) * sin(b * beta) * cos(th))
                worst = max(worst, abs(MU[0] - closed))
    good = worst < TOL
    n_checks += 1
    ok &= good
    print(f"    max |direct - closed form| over 5 knots x all arcs x 11 angles"
          f" = {worst:.2e}   {'PASS' if good else 'FAIL'}")

    # (2) the count ---------------------------------------------------------
    print("\n(2) N(3,n) by direct enumeration of admissible arcs, vs the closed form")
    print("     n   n%6   arcs  admissible  formula   l1=1 / l1=2 sheets")
    printed = {2: 1, 4: 3, 5: 4, 7: 4, 8: 5, 10: 7, 11: 8, 13: 8}
    for n in range(2, 26):
        if gcd(3, n) != 1:
            continue
        tc = traceless_characters(n)
        N, F = len(tc), N_formula(n)
        s1 = sum(1 for l1, _, _ in tc if l1 == 1)
        s2 = N - s1
        agree = (N == F)
        # Remark 4.4's symmetry (l1,l2) -> (3-l1, n-l2) preserves (*) and swaps
        # the sheets, but it preserves ADMISSIBILITY only for n odd (for n even
        # it breaks the parity match l1 = l2 mod 2).  So balance is an n-odd claim.
        sheets_ok = (s1 == s2) if n % 2 else True
        n_checks += 1 + (n % 2)
        ok &= agree and sheets_ok
        flag = "PASS" if agree else "FAIL"
        extra = ""
        if n in printed:
            match = (N == printed[n])
            n_checks += 1
            ok &= match
            extra = f"   paper says {printed[n]} {'ok' if match else 'MISMATCH'}"
        print(f"    {n:2d}    {n%6}    {len(arcs(n)):3d}      {N:3d}       {F:3d}"
              f"   {flag}   {s1} / {s2}{extra}")

    # (3) every admissible arc really is a traceless irreducible rep --------
    print("\n(3) Each admissible arc rebuilt as quaternions: rho(x)^3 = rho(y)^n,"
          " tr rho(mu) = 0")
    worst_tr, worst_z, sheet_minus, sheet_plus = 0.0, 0.0, 0, 0
    for n in range(2, 26):
        if gcd(3, n) != 1:
            continue
        a, b = meridian_exponents(n)
        for l1, l2, c in traceless_characters(n):
            th = math.acos(max(-1.0, min(1.0, c)))
            alpha, beta = pi * l1 / 3.0, pi * l2 / n
            X = qexp(alpha, (0.0, 0.0, 1.0))
            Y = qexp(beta, (sin(th), 0.0, cos(th)))
            X3, Yn = qpow(X, 3), qpow(Y, n)
            worst_z = max(worst_z, max(abs(X3[i] - Yn[i]) for i in range(4)))
            sgn = (-1) ** l1
            worst_z = max(worst_z, abs(X3[0] - sgn),
                          max(abs(X3[i]) for i in range(1, 4)))
            MU = qmul(qpow(X, a), qpow(Y, b))
            worst_tr = max(worst_tr, abs(2.0 * MU[0]))
            sheet_minus += (sgn == -1)
            sheet_plus += (sgn == 1)
    good_tr, good_z = worst_tr < 1e-9, worst_z < 1e-9
    n_checks += 2
    ok &= good_tr and good_z
    print(f"    max |tr rho(mu)|                  = {worst_tr:.2e}"
          f"   {'PASS' if good_tr else 'FAIL'}")
    print(f"    max |rho(x)^3 - rho(y)^n| and     = {worst_z:.2e}"
          f"   {'PASS' if good_z else 'FAIL'}")
    print(f"      |rho(z) - (-1)^l1|")
    print(f"    central sign: rho(z) = -1 on {sheet_minus} characters (l1 = 1),"
          f"  +1 on {sheet_plus} (l1 = 2)")
    print( "      -> Sec. 4.2's 'rho(x)^3 = rho(y)^n = -1' is the l1 = 1 sheet;"
           " the l1 = 2 sheet has +1.")

    # (4) Theorem B: the dihedral dichotomy ---------------------------------
    print("\n(4) Theorem B: irreducible traceless DIHEDRAL characters = (det - 1)/2")
    for n in range(2, 26):
        if gcd(3, n) != 1:
            continue
        # rho(y) traceless <=> beta = pi/2 <=> l2 = n/2; rho(x) never traceless
        dihedral = sum(1 for l1, l2, _ in traceless_characters(n) if 2 * l2 == n)
        x_traceless = sum(1 for l1 in (1, 2) if abs(cos(pi * l1 / 3.0)) < 1e-14)
        det = 3 if n % 2 == 0 else 1
        expect = (det - 1) // 2
        good = (dihedral == expect) and (x_traceless == 0)
        n_checks += 1
        ok &= good
        if n <= 14 or not good:
            print(f"    n = {n:2d}: det = {det}, dihedral = {dihedral},"
                  f" expected {expect}   {'PASS' if good else 'FAIL'}")
    print("    (all n <= 25 checked; n odd => 0 dihedral => every irreducible"
          " traceless character is non-dihedral)")

    # (5) N(3,n) = 2a for n odd, against fs_gradings.py ---------------------
    print("\n(5) N(3,n) = 2a for all odd n <= 43, a = -sigma/4 = #(spherical"
          " triangles), vs fs_gradings.py")
    try:
        from fs_gradings import valid
    except ImportError:
        print("    fs_gradings.py not alongside -- SKIPPED")
        valid = None
    if valid is not None:
        line = []
        for n in range(5, 44, 2):
            if gcd(3, n) != 1:
                continue
            a_count = len([l3 for l3 in range(2, n, 2) if valid(n, l3)])
            N = len(traceless_characters(n))
            good = (N == 2 * a_count)
            n_checks += 1
            ok &= good
            line.append(f"n={n}: N={N}, 2a={2*a_count} {'ok' if good else 'FAIL'}")
        for i in range(0, len(line), 3):
            print("    " + ";  ".join(line[i:i + 3]))
        print("    => sigma(T(3,n)) = -2 N(3,n) for all odd n <= 43")

    print(f"\n{'ALL PASS' if ok else 'FAILURES PRESENT'}  ({n_checks} checks)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
