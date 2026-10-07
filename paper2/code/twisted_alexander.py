#!/usr/bin/env python3
r"""Section 5.2-5.3 and Appendix A.1(iii): Alexander polynomials, determinants and signatures of the twisted torus
knots T(3, n; 2, m), the closures of (sigma_1 sigma_2)^n sigma_1^(2m).

Reduced Burau representation of B_3: psi(sigma_1) = [[-t, 1], [0, 1]], psi(sigma_2) = [[1, 0], [t, -t]], and
Delta of a closed 3-braid knot = det(I - psi(beta)) / (1 + t + t^2) up to a unit.  The script checks:

  (1) the matrix facts in the proof of Proposition 5.4: the braid relation, P = psi(sigma_1 sigma_2) has trace -t and
      determinant t^2 and P^3 = t^3 I, psi(sigma_1)^j = [[(-t)^j, (1 - (-t)^j)/(1 + t)], [0, 1]], and the trace and
      determinant of X = t^(3 nu) P psi(sigma_1)^(2M);
  (2) formula (5.1) for F_{nu,M} and the three telescoping identities of the proof, as identities of rational
      functions in independent variables t, W = t^(3 nu), T = t^(2M), together with the case n = 3 nu + 2;
  (3) Lemma 5.3 in the Burau representation, for nu, m <= 10 (the lemma is proved by hand);
  (4) Proposition 5.4 exactly: for every n <= 80 prime to 3 and 0 <= m <= 20 (1092 knots) the Burau Alexander
      polynomial equals the closed form, its coefficients lie in {-1, 0, 1}, its degree is 2(n + m - 1), and
      ell = 4 nu + 2M + 1; for n = 1, 2 (mod 6), n >= 7, this equals |sigma| + 1 with sigma = -(4 nu + 2M)
      (Proposition 5.5); the example T(3,7;2,2) of Section 5.2;
  (5) Proposition 5.5: the matrices at t = -1 in its proof, symbolically in m, and det T(3,n;2,m) = 2m + 1
      (n = 1 mod 6) or 2m + 3 (n = 2 mod 6) for 7 <= n <= 50 and 0 <= m <= 20; |sigma(T(3,n))| = 2(n-1) - 4 floor(n/6)
      equals 4 nu or 4 nu + 2;
  (6) the case n = 5 (Theorem 4.4): the closure of (sigma_1 sigma_2)^5 sigma_1^(2m) has the Alexander polynomial (3.2)
      of P(-2, 3, 5 + 2m), with ell = q + 2 and Delta(-1) = 6 - q, for 0 <= m <= 23;
  (7) the inputs of Lemma 5.6 for 2 <= nu <= 20, 2 <= M <= 20: the coefficients of t^2 and t^5 vanish, more than three
      consecutive coefficients are nonzero, and det = |2 + (-1)^nu (2M - 1)| < ell; the coefficient of t^5 in (3.2) is 1.

All computations are exact.
"""
import sys
import time

import sympy as sp

t, W, T, X, m_ = sp.symbols('t W T X m')
s1 = sp.Matrix([[-t, 1], [0, 1]])
s2 = sp.Matrix([[1, 0], [t, -t]])
P = s1 * s2
checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


def zero(expr):
    return sp.simplify(sp.together(expr)) == 0


def mzero(M):
    return sp.simplify(M) == sp.zeros(*M.shape)


def normalize(coeffs):
    c = list(coeffs)
    while c and c[0] == 0:
        c = c[1:]
    while c and c[-1] == 0:
        c = c[:-1]
    if c and c[0] < 0:
        c = [-x for x in c]
    return tuple(c)


def closed_form(nu, M):
    """Proposition 5.4: the coefficients of L + L_mid + L^* (lowest degree first)."""
    d = {}
    for i in range(nu):
        d[3 * i] = d.get(3 * i, 0) + 1
        d[3 * i + 1] = d.get(3 * i + 1, 0) - 1
    for i in range(2 * M + 1):
        d[3 * nu + i] = d.get(3 * nu + i, 0) + (-1) ** i
    top = 6 * nu + 2 * M
    for i in range(nu):
        d[top - 3 * i] = d.get(top - 3 * i, 0) + 1
        d[top - 3 * i - 1] = d.get(top - 3 * i - 1, 0) - 1
    return [d.get(e, 0) for e in range(max(d) + 1)]


def burau_delta(n, m):
    base = sp.eye(2) if n % 3 == 0 else (P if n % 3 == 1 else P ** 2)
    S = sp.Matrix([[t ** (2 * m), sp.cancel((1 - t ** (2 * m)) / (1 + t))], [0, 1]])
    Bm = (t ** (3 * (n // 3)) * base * S).applyfunc(sp.expand)
    det = sp.Poly(sp.expand((sp.eye(2) - Bm).det()), t)
    q, r = sp.div(det, sp.Poly(1 + t + t ** 2, t))
    assert r.is_zero, (n, m)
    return normalize([int(c) for c in q.all_coeffs()[::-1]])


def nuM(n, m):
    return ((n - 1) // 3, m) if n % 3 == 1 else ((n - 2) // 3, m + 1)


def main():
    t0 = time.time()
    print('(1) Burau matrices (proof of Proposition 5.4)')
    check('braid relation psi(s1)psi(s2)psi(s1) = psi(s2)psi(s1)psi(s2)', mzero(s1 * s2 * s1 - s2 * s1 * s2))
    check('P = psi(s1 s2) = [[0,-t],[t,-t]], tr P = -t, det P = t^2, P^3 = t^3 I',
          mzero(P - sp.Matrix([[0, -t], [t, -t]])) and zero(P.trace() + t) and zero(P.det() - t ** 2)
          and mzero(P ** 3 - t ** 3 * sp.eye(2)))
    Sj = sp.Matrix([[X, (1 - X) / (1 + t)], [0, 1]])          # X stands for (-t)^j
    check('psi(s1)^j = [[(-t)^j, (1-(-t)^j)/(1+t)],[0,1]] (induction step, and j <= 40 directly)',
          mzero(Sj * s1 - sp.Matrix([[-t * X, (1 + t * X) / (1 + t)], [0, 1]]))
          and all(mzero(s1 ** j - sp.Matrix([[(-t) ** j, sp.cancel((1 - (-t) ** j) / (1 + t))], [0, 1]]))
                  for j in range(41)))
    ST = sp.Matrix([[T, (1 - T) / (1 + t)], [0, 1]])          # psi(s1)^(2M), T = t^(2M)
    Xm = W * P * ST                                           # W = t^(3 nu)
    check('tr X = -t^(3nu+1) (t + T)/(1 + t) and det X = t^(6nu+2) T',
          zero(Xm.trace() + t * W * (t + T) / (1 + t)) and zero(Xm.det() - t ** 2 * W ** 2 * T))

    print('(2) formula (5.1) and the telescoping identities (W = t^(3nu), T = t^(2M) independent)')
    F = 1 + t ** 2 * W * (1 + T / t) / (1 + t) + t ** 2 * W ** 2 * T
    check('det(I - X) = 1 + t^(3nu+2)(1 + t^(2M-1))/(1+t) + t^(6nu+2M+2)', zero((sp.eye(2) - Xm).det() - F))
    check('M = 0: the middle term is t^(3nu+1)', zero(F.subs(T, 1) - (1 + t * W + t ** 2 * W ** 2)))
    F2 = (sp.eye(2) - W * P ** 2 * ST).det()
    check('n = 3nu+2: det(I - t^(3nu) P^2 psi(s1)^(2M)) = F with T replaced by t^2 T (M -> M + 1)',
          zero(F2 - F.subs(T, t ** 2 * T)))
    L = (1 - t) * (1 - W) / (1 - t ** 3)
    Lstar = W ** 2 * T * L.subs({t: 1 / t, W: 1 / W})
    Lmid = W * (1 + t * T) / (1 + t)
    check('(1+t+t^2) L = 1 - t^(3nu)', zero((1 + t + t ** 2) * L - (1 - W)))
    check('(1+t+t^2) L^* = t^(6nu+2M+2) - t^(3nu+2M+2)', zero((1 + t + t ** 2) * Lstar - (t ** 2 * W ** 2 * T - t ** 2 * W * T)))
    check('(1+t+t^2)(1+t^(2M+1)) - (1+t)(1+t^(2M+2)) = t^2 + t^(2M+1)',
          zero((1 + t + t ** 2) * (1 + t * T) - (1 + t) * (1 + t ** 2 * T) - (t ** 2 + t * T)))
    check('(1+t+t^2) L_mid = t^(3nu) + t^(3nu+2M+2) + t^(3nu+2)(1+t^(2M-1))/(1+t)',
          zero((1 + t + t ** 2) * Lmid - (W + t ** 2 * W * T + t ** 2 * W * (1 + T / t) / (1 + t))))
    check('(1+t+t^2)(L + L_mid + L^*) = F', zero((1 + t + t ** 2) * (L + Lmid + Lstar) - F))

    print('(3) Lemma 5.3 in the Burau representation, nu, m <= 10')
    check('psi(s1^-1 (s1s2)^(3nu+2) s1^(2m) s1) = psi((s1s2)^(3nu+1) s1^(2m+2))',
          all(mzero(s1.inv() * P ** (3 * nu + 2) * s1 ** (2 * m) * s1 - P ** (3 * nu + 1) * s1 ** (2 * m + 2))
              for nu in range(11) for m in range(11)))

    print('(4) Proposition 5.4 exactly: n <= 80 prime to 3, m <= 20')
    t1 = time.time()
    count, ok_cf, ok_coef, ok_deg, ok_norm, ok_sig = 0, True, True, True, True, True
    for n in [n for n in range(4, 81) if n % 3]:
        for m in range(21):
            nu, M = nuM(n, m)
            D = burau_delta(n, m)
            ok_cf &= D == normalize(closed_form(nu, M))
            ok_coef &= set(D) <= {-1, 0, 1}
            ok_deg &= len(D) - 1 == 2 * (n + m - 1)
            ell = sum(abs(c) for c in D)
            ok_norm &= ell == 4 * nu + 2 * M + 1
            if n % 6 in (1, 2) and n >= 7:
                sig = -(2 * (n - 1) - 4 * (n // 6)) - 2 * m          # sigma(T(3,n)) - 2m, Proposition 5.5
                ok_sig &= sig == -(4 * nu + 2 * M) and ell == abs(sig) + 1
            count += 1
    print(f'      {count} knots ({time.time() - t1:.0f} s)')
    check(f'{count} knots (paper: 1092): Burau Alexander polynomial = closed form of Proposition 5.4', count == 1092 and ok_cf)
    check('coefficients in {-1, 0, 1}; degree 6nu + 2M = 2(n + m - 1)', ok_coef and ok_deg)
    check('ell = 4 nu + 2M + 1; for n = 1, 2 (mod 6), n >= 7: sigma(T(3,n)) - 2m = -(4nu + 2M) and ell = |sigma| + 1',
          ok_norm and ok_sig)
    ex = burau_delta(7, 2)
    paper_ex = (1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 1, 0, -1, 1, 0, -1, 1)
    check(f'T(3,7;2,2): coefficients of t^0..t^16 = {ex}, ell = 13 = |sigma| + 1', ex == paper_ex and sum(map(abs, ex)) == 13)

    print('(5) Proposition 5.5: determinants and signatures')
    s1m, s2m = s1.subs(t, -1), s2.subs(t, -1)
    Q = s1m * s2m
    check('Q = psi(s1 s2)|_{t=-1} = [[0,1],[-1,1]] has order 6',
          Q == sp.Matrix([[0, 1], [-1, 1]]) and Q ** 6 == sp.eye(2) and all(Q ** k != sp.eye(2) for k in range(1, 6)))
    S2m = sp.Matrix([[1, 2 * m_], [0, 1]])                    # psi(s1)^(2m) at t = -1
    check('psi(s1)^(2m)|_{t=-1} = [[1, 2m],[0,1]] (m <= 20)', all(s1m ** (2 * k) == S2m.subs(m_, k) for k in range(21)))
    M1, M2 = Q * S2m, Q ** 2 * S2m
    check('n = 1 (mod 6): Q^n psi(s1)^(2m) = [[0,1],[-1,1-2m]], det(I - .) = 2m + 1',
          sp.simplify(M1 - sp.Matrix([[0, 1], [-1, 1 - 2 * m_]])) == sp.zeros(2) and sp.expand((sp.eye(2) - M1).det() - (2 * m_ + 1)) == 0)
    check('n = 2 (mod 6): Q^n psi(s1)^(2m) = [[-1,1-2m],[-1,-2m]], det(I - .) = 2m + 3',
          sp.simplify(M2 - sp.Matrix([[-1, 1 - 2 * m_], [-1, -2 * m_]])) == sp.zeros(2) and sp.expand((sp.eye(2) - M2).det() - (2 * m_ + 3)) == 0)
    ok = True
    for n in [n for n in range(7, 51) if n % 6 in (1, 2)]:
        for m in range(21):
            d = abs((sp.eye(2) - Q ** n * s1m ** (2 * m)).det())
            ok &= d == (2 * m + 1 if n % 6 == 1 else 2 * m + 3)
    check('det T(3,n;2,m) = 2m+1 (n = 1) or 2m+3 (n = 2 mod 6) for 7 <= n <= 50, m <= 20', ok)
    check('|sigma(T(3,n))| = 2(n-1) - 4 floor(n/6) = 4 nu (n = 6j+1), 4 nu + 2 (n = 6j+2), n <= 200',
          all(2 * (n - 1) - 4 * (n // 6) == (4 * (n // 3) if n % 6 == 1 else 4 * (n // 3) + 2)
              for n in range(7, 201) if n % 6 in (1, 2)))

    print('(6) n = 5: P(-2,3,5+2m) (Lemma 3.1, Theorem 4.4), m <= 23')
    ok, t5 = True, True
    for m in range(24):
        q = 5 + 2 * m
        D = burau_delta(5, m)
        claim = normalize([1, -1, 0] + [(-1) ** (j + 1) for j in range(3, q + 1)] + [0, -1, 1])
        ok &= D == claim and sum(map(abs, D)) == q + 2 and abs(sum(c * (-1) ** i for i, c in enumerate(D))) == abs(6 - q)
        t5 &= claim[5] == 1
    check('closure of (s1s2)^5 s1^(2m) has Delta = (3.2) with q = 5 + 2m, ell = q + 2, |Delta(-1)| = |6 - q|', ok)

    print('(7) inputs of Lemma 5.6, 2 <= nu, M <= 20')
    ok = True
    for nu in range(2, 21):
        for M in range(2, 21):
            D = closed_form(nu, M)
            ell = sum(map(abs, D))
            run, best = 0, 0
            for c in D:
                run = run + 1 if c else 0
                best = max(best, run)
            det = abs(sum(c * (-1) ** i for i, c in enumerate(D)))
            ok &= (D[2] == 0 and D[5] == 0 and best >= 2 * M + 1 > 3 and det == abs(2 + (-1) ** nu * (2 * M - 1)) and det < ell)
    check('coefficients of t^2 and t^5 vanish; 2M + 1 >= 5 consecutive nonzero coefficients; det = |2 + (-1)^nu (2M-1)| < ell', ok)
    check('torus knots T(3,b) = K_{nu,0} or K_{nu,1}: at most three consecutive nonzero coefficients (nu <= 30)',
          all(max(len(r) for r in ''.join('x' if c else ' ' for c in closed_form(nu, M)).split()) <= 3
              for nu in range(1, 31) for M in (0, 1)))
    check('coefficient of t^5 in (3.2) is 1 for 5 <= q <= 51', t5)

    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
