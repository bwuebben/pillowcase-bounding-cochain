#!/usr/bin/env python3
r"""Section 4 and Appendix A.3(ii): the combinatorics of the sheared polygons, in exact rational arithmetic, for every
odd q with 7 <= q <= 401 (or up to the bound given on the command line).

Coordinates are in units of pi.  The lift of the unperturbed curve of the T(3,5) tangle (Corollary 4.8) is the polygon

    O = (0,0) -B_0-> J_1 = (1/6, 0) -H-> J_2 = (5/6, 4) -B_1^-1-> J_3 = (1/6, 4) -H'-> J_4 = (5/6, 8) -B_2-> E = (1, 8),

with E' = (1, 4) when B_2 follows H (split case), and the remnant circle B_1^-1 . H' is (0,4)-periodic.  The shear is
S_{-m}(gamma, theta) = (gamma, theta + (q-5) gamma), phi = theta - gamma, and Delta_k = {phi = 2k}.  For each q the script
checks:

  * the values (4.5) of phi at the sheared vertices, and that phi_1, ..., phi_4 lie at distance >= 1/6 from 2Z;
  * Remark 4.12: the distance from phi_3 to the next even integer above it is 1/6 exactly when q = 5 (mod 12), and the
    distance from phi_2 to the next even integer below it is 1/6 exactly when q = 11 (mod 12);
  * Proposition 4.13 (joined case): the lines crossed, three times (up, down, up) for k = 2 + j with
    (q-6)/12 < j < 5(q-6)/12 and once (up) for the other k in {1, ..., (q+1)/2}; N_q of them are crossed three times,
    and the arc carries q + 2 + 4 N_q generators; (q+2) - phi_2 = (q+18)/6;
  * Proposition 4.13 (split case): the arc crosses Delta_1, ..., Delta_{(q-3)/2} once and carries q - 2 generators;
    one period of the remnant circle crosses N_q lines downward and N_q + 2 upward, so the circle carries 4 + 4 N_q;
  * Proposition 4.13, last paragraph: every crossing abscissa is 2i/(q-6) or (2i+1)/q, never 1/2;
  * Proposition 4.17: on a line crossed three times the abscissae are gamma_1 = (2k+1)/q, gamma_2 = 2(k-2)/(q-6),
    gamma_3 = (2k-3)/q, the two differences are the displayed fractions and are >= 1/(q(q-6)), and
    5/6 > gamma_1 > gamma_2 > gamma_3 > 1/6;
  * Proposition 4.15 and Theorem 4.18(3): the degrees 2k (up) and 2k+1 (down) of x^+, and 2k-1, 2k of x^-; the graded
    homology of the strand complexes of Proposition 4.17 equals (1 + N_3, N_1, N_1, N_3) in both resolutions, and the
    rank of the differential is 2 N_q;
  * Remark 4.19: q + 2 + 4 N_q = q - 6 + 4 (floor((5q+6)/12) - floor((q+6)/12)).

These are exact computations.  They check, for each q in the range, the combinatorial steps of the hand proofs of
Propositions 4.13 and 4.17 and Theorem 4.18.
"""
import sys
import time
from fractions import Fraction as F
from math import floor

checks = []
failures = []


def check(desc, ok, quiet=False):
    ok = bool(ok)
    checks.append(ok)
    if not ok:
        failures.append(desc)
    if not quiet or not ok:
        print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def Nq(q):
    return floor(F(5 * (q - 6), 12)) - floor(F(q - 6, 12))


def crossings(P, Q, q):
    """Crossings of the sheared segment P -> Q with the lines Delta_k, k >= 1: (k, gamma, 'up' or 'down')."""
    (g0, t0), (g1, t1) = P, Q
    f0, f1 = t0 + (q - 6) * g0, t1 + (q - 6) * g1          # phi o S_{-m} = theta + (q-6) gamma
    out = []
    lo, hi = min(f0, f1), max(f0, f1)
    k = floor(lo / 2) + 1
    while 2 * k < hi:
        if 2 * k > lo and k >= 1:
            lam = (2 * k - f0) / (f1 - f0)
            out.append((k, g0 + lam * (g1 - g0), 'up' if f1 > f0 else 'down'))
        k += 1
    return out


def graded(gens_by_line):
    """Graded homology (degrees mod 4) of the strand complexes of Proposition 4.17, plus r_+ in degree 0."""
    H = [1, 0, 0, 0]
    rank = 0
    for k, dirs in gens_by_line.items():
        # on S^k_+ (x^+) the degrees are 2k (up), 2k+1 (down); on S^k_- one less
        for sign_shift in (0, -1):
            degs = [((2 * k if d == 'up' else 2 * k + 1) + sign_shift) % 4 for d in dirs]
            if dirs == ['up']:
                H[degs[0]] += 1
            elif dirs == ['up', 'down', 'up']:
                # d x_{P2} = x_{P1} + x_{P3}: rank one, homology in the degree of the up-crossings
                assert degs[0] == degs[2] == (degs[1] - 1) % 4
                H[degs[0]] += 1
                rank += 1
            else:
                raise AssertionError(dirs)
    return H, rank


def one_q(q, verbose):
    V = [(F(0), F(0)), (F(1, 6), F(0)), (F(5, 6), F(4)), (F(1, 6), F(4)), (F(5, 6), F(8)), (F(1), F(8))]
    Ep = (F(1), F(4))
    phi = [t + (q - 6) * g for (g, t) in V]
    phiEp = Ep[1] + (q - 6) * Ep[0]
    claim = [F(0), F(q - 6, 6), 4 + F(5 * (q - 6), 6), 4 + F(q - 6, 6), 8 + F(5 * (q - 6), 6), F(q + 2)]
    check(f'q={q}: values (4.5) of phi at O, J_1, ..., J_4, E, E\'', phi == claim and phiEp == q - 2, quiet=not verbose)
    for i in (1, 2, 3, 4):
        r = phi[i] % 2
        check(f'q={q}: phi_{i} at distance >= 1/6 from 2Z', min(r, 2 - r) >= F(1, 6), quiet=True)
    up3 = (-phi[3]) % 2                                  # distance from phi_3 up to the next even integer
    dn2 = phi[2] % 2                                     # distance from phi_2 down to the next even integer
    check(f'q={q}: Remark 4.12 margins', up3 >= F(1, 6) and dn2 >= F(1, 6)
          and (up3 == F(1, 6)) == (q % 12 == 5) and (dn2 == F(1, 6)) == (q % 12 == 11), quiet=not verbose)
    # joined case
    cr = []
    for i in range(5):
        cr += [(i,) + c for c in crossings(V[i], V[i + 1], q)]
    lines = {}
    for (i, k, g, d) in cr:
        lines.setdefault(k, []).append((i, g, d))
    triple = sorted(k for k in lines if len(lines[k]) == 3)
    pred_triple = sorted(2 + j for j in range(-10, q) if F(q - 6, 12) < j < F(5 * (q - 6), 12))
    check(f'q={q}: joined arc crosses exactly Delta_1, ..., Delta_{(q + 1) // 2}, each once or three times',
          sorted(lines) == list(range(1, (q + 1) // 2 + 1)) and all(len(v) in (1, 3) for v in lines.values()),
          quiet=not verbose)
    check(f'q={q}: triple lines = {{2 + j : (q-6)/12 < j < 5(q-6)/12}}, N_q = {Nq(q)} of them',
          triple == pred_triple and len(triple) == Nq(q), quiet=not verbose)
    check(f'q={q}: joined arc carries 1 + 2 #crossings = q + 2 + 4 N_q = {q + 2 + 4 * Nq(q)} generators',
          1 + 2 * len(cr) == q + 2 + 4 * Nq(q), quiet=not verbose)
    check(f'q={q}: (q + 2) - phi_2 = (q + 18)/6', (q + 2) - phi[2] == F(q + 18, 6), quiet=True)
    for k in triple:
        (i1, g1, d1), (i2, g2, d2), (i3, g3, d3) = lines[k]
        x = 2 * (k - 2)
        ok = ((i1, i2, i3) == (1, 2, 3) and (d1, d2, d3) == ('up', 'down', 'up')
              and g1 == F(2 * k + 1, q) and g2 == F(2 * (k - 2), q - 6) and g3 == F(2 * k - 3, q)
              and g2 - g3 == 6 * (x - F(q - 6, 6)) / (q * (q - 6)) and g1 - g2 == 6 * (F(5 * (q - 6), 6) - x) / (q * (q - 6))
              and min(g2 - g3, g1 - g2) >= F(1, q * (q - 6)) and F(5, 6) > g1 > g2 > g3 > F(1, 6))
        check(f'q={q}, k={k}: Proposition 4.17 abscissae, differences and order', ok, quiet=True)
    absc = [g for (_, _, g, _) in cr]
    check(f'q={q}: every abscissa is 2i/(q-6) or (2i+1)/q, never 1/2',
          all(((g * (q - 6) / 2).denominator == 1) or ((g * q - 1) / 2).denominator == 1 for g in absc)
          and all(g != F(1, 2) for g in absc), quiet=True)
    gens_by_line = {k: [d for (_, _, d) in v] for k, v in lines.items()}
    H, rank = graded(gens_by_line)
    pred = [1 + (q + 1) // 4, -(-(q + 1) // 4), -(-(q + 1) // 4), (q + 1) // 4]
    check(f'q={q}: joined homology {tuple(H)} = (1 + N_3, N_1, N_1, N_3), rank of d = 2 N_q = {2 * Nq(q)}',
          H == pred and rank == 2 * Nq(q), quiet=not verbose)
    # split case: the arc O -> J_1 -> J_2 -> E'
    arc_split = crossings(V[0], V[1], q) + crossings(V[1], V[2], q) + crossings(V[2], Ep, q)
    ks = sorted(k for (k, _, _) in arc_split)
    check(f'q={q}: split arc crosses Delta_1, ..., Delta_{(q - 3) // 2} once (up), q - 2 = {q - 2} generators',
          ks == list(range(1, (q - 3) // 2 + 1)) and all(d == 'up' for (_, _, d) in arc_split)
          and 1 + 2 * len(arc_split) == q - 2, quiet=not verbose)
    dec = crossings(V[2], V[3], q)                       # one period: B_1^-1 then H'
    inc = crossings(V[3], (V[4][0], V[4][1]), q)
    check(f'q={q}: one period of the remnant circle crosses N_q lines down and N_q + 2 up; 4 + 4 N_q = '
          f'{4 + 4 * Nq(q)} generators', len(dec) == Nq(q) and len(inc) == Nq(q) + 2
          and all(d == 'down' for (_, _, d) in dec) and all(d == 'up' for (_, _, d) in inc)
          and 2 * (len(dec) + len(inc)) == 4 + 4 * Nq(q), quiet=not verbose)
    Hs = [1, 0, 0, 0]
    for k in ks:
        Hs[(2 * k) % 4] += 1
        Hs[(2 * k - 1) % 4] += 1
    check(f'q={q}: split homology = arc + (1,1,1,1) = (1 + N_3, N_1, N_1, N_3)',
          [Hs[i] + 1 for i in range(4)] == pred, quiet=not verbose)
    check(f'q={q}: Remark 4.19 identity',
          q + 2 + 4 * Nq(q) == q - 6 + 4 * (floor(F(5 * q + 6, 12)) - floor(F(q + 6, 12))), quiet=True)


def main():
    qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 401
    t0 = time.time()
    print('Values of Section 4.5 and Theorem 4.18')
    table = {7: 0, 9: 1, 11: 2, 13: 2, 15: 3, 17: 4, 19: 4}
    for q, n in table.items():
        check(f'N_{q} = {n}', Nq(q) == n)
    check('q = 7: 9 = 1 + |sigma(K_7)| generators (sigma = -(q+1))', 7 + 2 + 4 * Nq(7) == 1 + 8)
    print(f'Exact checks for odd q in [7, {qmax}] (printed in full for q <= 13)')
    for q in range(7, qmax + 1, 2):
        one_q(q, verbose=(q <= 13))
    print('  q    N_q  generators  rank d  homology (deg 0..3)')
    for q in (7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 33):
        n = Nq(q)
        print(f'  {q:3d}  {n:3d}  {q + 2 + 4 * n:10d}  {2 * n:6d}  '
              f'({1 + (q + 1) // 4}, {-(-(q + 1) // 4)}, {-(-(q + 1) // 4)}, {(q + 1) // 4})')
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    for f in failures[:20]:
        print('  failed:', f)
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
