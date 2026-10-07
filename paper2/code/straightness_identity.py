#!/usr/bin/env python3
r"""Lemma 4.7 and Appendix A.1(i): the image of the central circle C_0 of the T(3,5) tangle is straight.

The tangle is that of Hedden, Herald and Kirk for T(3,5) with (r, s) = (2, -1).  Its group is free on A, B, and at
eps = 0 the boundary meridians are (4.1)

    a = A^2 B^3,   b = B^-2 A^-1,   c = B^-2 a B^2.

Lemma 4.7 states that w = (ba)^-3 c (ba)^3 a satisfies rho(w) = 1 on C_0, so that theta = 6 gamma - pi on Pi(C_0).
This script checks:

  (1) every displayed identity in Steps 1-7 of the proof of Lemma 4.7, symbolically in the indeterminates v, zeta
      (quaternions written p = p_0 + vec p; beta = rho(b), N = rho(B) = cos v + sin v Q, zeta = beta . Q,
      U = cos v beta + sin v (beta x Q), E = 4 cos^2 v + 4 zeta^2 sin^2 v - 3, Q' = 2 zeta beta - Q, xi = beta N^-1);
  (2) the identities stated in Appendix A.1(i) without the assumption E = 0:
      vec rho(g) . vec rho(a) = (1/2) zeta sin 2v E^2  and  tr(rho(g) rho(a)) = -2 zeta sin 2v E^2;
  (3) the Fricke form: with x = tr A, y = tr B, z = tr AB one has tr b = yz - x, and on {x = yz}
      E = y^2 + z^2 - 3, tr a = -yE, tr g = -zE, tr(ga) = yz E^2, while tr w - 2 lies in the ideal (x - yz, E);
      for the word (ba)^-2 c (ba)^2 a the same reduction leaves -y^2 - 1;
  (4) rho(w) = 1 in exact quaternion arithmetic at three points of C_0 with algebraic coordinates;
  (5) at 50 significant digits, on the parametrization of C_0 in Lemma 3.3 and Definition 3.4 of the companion paper
      on T(3,n), |theta - 6 gamma + pi| <= 1.4e-47 modulo 2 pi, and gamma is strictly monotone along each half.

Items (1)-(4) are exact computations; (5) is a high-precision evaluation.  Lemma 4.7 is proved by hand in the paper;
these computations check the displayed steps.
"""
import sys
import time

import mpmath as mp
import sympy as sp

checks = []


def check(desc, ok):
    ok = bool(ok)
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {desc}", flush=True)


# ---------------------------------------------------------------------------------------------- quaternions
def qm(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return [w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2]


def conj(p):
    return [p[0], -p[1], -p[2], -p[3]]


def add(p, q, s=1):
    return [p[i] + s * q[i] for i in range(4)]


def scal(c, p):
    return [c * x for x in p]


def is_zero(expr):
    return sp.simplify(sp.expand(sp.expand_trig(expr))) == 0


def q_is_zero(p):
    return all(is_zero(x) for x in p)


def vec(p):
    return sp.Matrix(p[1:])


# ---------------------------------------------------------------------------------------------- free group
INV = {'A': 'a', 'a': 'A', 'B': 'b', 'b': 'B'}       # capital letter = generator, lower case = inverse


def inv(w):
    return ''.join(INV[ch] for ch in reversed(w))


def red(w):
    st = []
    for ch in w:
        if st and st[-1] == INV[ch]:
            st.pop()
        else:
            st.append(ch)
    return ''.join(st)


W_A, W_B = 'AABBB', 'bba'                          # a = A^2 B^3, b = B^-2 A^-1
W_C = 'bb' + W_A + 'BB'                            # c = B^-2 a B^2
W_BA = W_B + W_A
W_G = 'AB' * 3 + 'BB'                              # g = (AB)^3 B^2
W_W = inv(W_BA) * 3 + W_C + W_BA * 3 + W_A        # w = (ba)^-3 c (ba)^3 a


def part1():
    print('(1) Steps 1-7 of the proof of Lemma 4.7, symbolically in (v, zeta)')
    v, zt, gam, th = sp.symbols('v zeta gamma theta', real=True)
    cv, sv = sp.cos(v), sp.sin(v)
    # Step 1 (free group): ba = B^-2 (AB) B^2, (ba)^3 = B^-2 g, w = g^-1 a g a
    check('Step 1: ba = B^-2 (AB) B^2', red(W_BA) == red('bb' + 'AB' + 'BB'))
    check('Step 1: (ba)^3 = B^-2 g', red(W_BA * 3) == red('bb' + W_G))
    check('Step 1: w = g^-1 a g a as reduced words', red(W_W) == red(inv(W_G) + W_A + W_G + W_A))
    # Step 2: generic data, up to a rotation: Q = k, beta a pure unit with beta . Q = zeta
    Q = [0, 0, 0, 1]
    beta = [0, sp.sqrt(1 - zt ** 2), 0, zt]
    N = [cv, 0, 0, sv]
    Ni = conj(N)
    bv, Qv = vec(beta), vec(Q)
    Uv = cv * bv + sv * bv.cross(Qv)
    E = 4 * cv ** 2 + 4 * zt ** 2 * sv ** 2 - 3
    rA = qm(conj(beta), qm(Ni, Ni))                  # rho(A) = rho(b)^-1 N^-2, from b = B^-2 A^-1
    check('Step 2: rho(A) = -beta N^-2', q_is_zero(add(rA, scal(-1, qm(beta, qm(Ni, Ni))), -1)))
    # Step 3
    Qp = add(scal(2 * zt, beta), Q, -1)
    check("Step 3: beta Q beta^-1 = Q' = 2 zeta beta - Q", q_is_zero(add(qm(qm(beta, Q), conj(beta)), Qp, -1)))
    eQp = [sp.cos(2 * v)] + [-sp.sin(2 * v) * x for x in Qp[1:]]          # exp(-2v Q')
    check("Step 3: beta N^-2 beta = -exp(-2v Q')", q_is_zero(add(qm(qm(beta, qm(Ni, Ni)), beta), eQp)))
    ra = qm(qm(rA, rA), qm(N, qm(N, N)))
    check("Step 3: rho(a) = rho(A)^2 N^3 = -exp(-2v Q') N", q_is_zero(add(ra, qm(eQp, N))))
    s_, t_, ang = sp.symbols('s t omega', real=True)
    P_ = [0, 1, 0, 0]
    R_ = [0, sp.cos(ang), sp.sin(ang), 0]
    eP = [sp.cos(s_)] + [sp.sin(s_) * x for x in P_[1:]]
    eR = [sp.cos(t_)] + [sp.sin(t_) * x for x in R_[1:]]
    check('Step 3: Re(e^{sP} e^{tR}) = cos s cos t - (P.R) sin s sin t (P, R pure units)',
          is_zero(qm(eP, eR)[0] - (sp.cos(s_) * sp.cos(t_) - sp.cos(ang) * sp.sin(s_) * sp.sin(t_))))
    check("Step 3: Q'.Q = 2 zeta^2 - 1", is_zero(vec(Qp).dot(Qv) - (2 * zt ** 2 - 1)))
    check('Step 3: Re rho(a) = -cos v E', is_zero(ra[0] + cv * E))
    check('Step 3: vec rho(a) = -sin 3v Q + 2 zeta sin 2v U',
          q_is_zero([0] + list(vec(ra) - (-sp.sin(3 * v) * Qv + 2 * zt * sp.sin(2 * v) * Uv))))
    # Step 5
    xi = qm(beta, Ni)
    rAB = qm(rA, N)
    check('Step 5: rho(AB) = -xi, xi = beta N^-1', q_is_zero(add(rAB, xi)))
    check('Step 5: Re xi = zeta sin v', is_zero(xi[0] - zt * sv))
    check('Step 5: xi^2 = 2 (Re xi) xi - 1', q_is_zero(add(qm(xi, xi), add(scal(2 * xi[0], xi), [1, 0, 0, 0], -1), -1)))
    xi3 = qm(xi, qm(xi, xi))
    check('Step 5: xi^3 = (4 zeta^2 sin^2 v - 1) xi - 2 zeta sin v',
          q_is_zero(add(xi3, add(scal(4 * zt ** 2 * sv ** 2 - 1, xi), [2 * zt * sv, 0, 0, 0], -1), -1)))
    check('Step 5: xi N^2 = beta N', q_is_zero(add(qm(xi, qm(N, N)), qm(beta, N), -1)))
    rg = qm(qm(qm(rAB, rAB), rAB), qm(N, N))
    check('Step 5: rho(g) = -xi^3 N^2', q_is_zero(add(rg, qm(xi3, qm(N, N)))))
    check('Step 5: rho(g) = -(4 zeta^2 sin^2 v - 1) beta N + 2 zeta sin v N^2',
          q_is_zero(add(rg, add(scal(-(4 * zt ** 2 * sv ** 2 - 1), qm(beta, N)), scal(2 * zt * sv, qm(N, N))), -1)))
    check('Step 5: Re rho(g) = zeta sin v E', is_zero(rg[0] - zt * sv * E))
    zeta_on_E = sp.sqrt(3 - 4 * cv ** 2) / (2 * sv)                     # a branch of E = 0
    check('Step 5: on E = 0, 4 zeta^2 sin^2 v - 1 = -2 cos 2v',
          is_zero((4 * zt ** 2 * sv ** 2 - 1 + 2 * sp.cos(2 * v)).subs(zt, zeta_on_E)))
    gvec_general = -(4 * zt ** 2 * sv ** 2 - 1) * Uv + 2 * zt * sv * sp.sin(2 * v) * Qv
    check('Step 5: vec rho(g) = -(4 zeta^2 sin^2 v - 1) U + 2 zeta sin v sin 2v Q (identically)',
          q_is_zero([0] + list(vec(rg) - gvec_general)))
    gvec_E = 2 * sp.cos(2 * v) * Uv + 2 * zt * sv * sp.sin(2 * v) * Qv
    check('Step 5: on E = 0, vec rho(g) = 2 cos 2v U + 2 zeta sin v sin 2v Q',
          q_is_zero([0] + [x.subs(zt, zeta_on_E) for x in (vec(rg) - gvec_E)]))
    # Step 6
    check('Step 6: U.U = 1 - zeta^2 sin^2 v', is_zero(Uv.dot(Uv) - (1 - zt ** 2 * sv ** 2)))
    check('Step 6: U.Q = zeta cos v', is_zero(Uv.dot(Qv) - zt * cv))
    check('Step 6: sin 3v = sin v (4 cos^2 v - 1)', is_zero(sp.sin(3 * v) - sv * (4 * cv ** 2 - 1)))
    avec = -sp.sin(3 * v) * Qv + 2 * zt * sp.sin(2 * v) * Uv
    bracket = 2 * zt ** 2 * sv ** 2 * sp.sin(2 * v) + 2 * sp.cos(2 * v) * sp.sin(2 * v) - sp.sin(3 * v) * cv
    check('Step 6: (2 cos 2v U + 2 zeta sin v sin 2v Q).(vec rho(a)) = 2 zeta (2 zeta^2 sin^2 v sin 2v '
          '+ 2 cos 2v sin 2v - sin 3v cos v)', is_zero(gvec_E.dot(avec) - 2 * zt * bracket))
    check('Step 6: 2 zeta (...) = zeta sin 2v E', is_zero(2 * zt * bracket - zt * sp.sin(2 * v) * E))
    # Step 7
    I = [0, 1, 0, 0]

    def ek(t):
        return [sp.cos(t), 0, 0, sp.sin(t)]
    ra_, rb_, rc_ = I, qm(ek(gam), I), qm(ek(th), I)
    rba = qm(rb_, ra_)
    check('Step 7: rho(ba) = -e^{gamma k}', q_is_zero(add(rba, ek(gam))))
    tt = sp.symbols('t', real=True)
    check('Step 7: i e^{tk} = e^{-tk} i', q_is_zero(add(qm(I, ek(tt)), qm(ek(-tt), I), -1)))
    for j in (-3, -2, -1, 1, 2, 3):
        P = [1, 0, 0, 0]
        base = rba if j > 0 else conj(rba)
        for _ in range(abs(j)):
            P = qm(P, base)
        word = qm(qm(qm(P, rc_), conj(P)), ra_)
        diff = add(word, scal(-1, ek(th + 2 * j * gam)), -1)
        check(f'Step 7: rho((ba)^{j} c (ba)^{-j} a) = -e^((theta + {2 * j} gamma) k)',
              all(sp.simplify(sp.expand(x.rewrite(sp.exp))) == 0 for x in diff))
    # (2) the general identities (Appendix A.1(i))
    print('(2) without the assumption E = 0 (Appendix A.1(i))')
    check('vec rho(g) . vec rho(a) = (1/2) zeta sin 2v E^2', is_zero(vec(rg).dot(vec(ra)) - zt * sp.sin(2 * v) * E ** 2 / 2))
    check('tr(rho(g) rho(a)) = -2 zeta sin 2v E^2', is_zero(2 * qm(rg, ra)[0] + 2 * zt * sp.sin(2 * v) * E ** 2))
    check('tr rho(B) = 2 cos v and tr rho(AB) = -2 zeta sin v, so y^2 + z^2 - 3 = E',
          is_zero(2 * N[0] - 2 * cv) and is_zero(2 * rAB[0] + 2 * zt * sv)
          and is_zero((2 * cv) ** 2 + (2 * zt * sv) ** 2 - 3 - E))


# ---------------------------------------------------------------------------------------------- Fricke form
def part3():
    print('(3) Fricke coordinates x = tr A, y = tr B, z = tr AB (Appendix A.1(i))')
    x, y, z, t = sp.symbols('x y z t')
    A = sp.Matrix([[x, -1], [1, 0]])
    B = sp.Matrix([[0, t], [-1 / t, y]])                 # tr B = y, tr AB = t + 1/t
    gen = {'A': A, 'a': A.inv(), 'B': B, 'b': B.inv()}

    def fricke(word):
        """The polynomial F with tr(word) = F(tr A, tr B, tr AB)."""
        M = sp.eye(2)
        for ch in word:
            M = (M * gen[ch]).applyfunc(sp.expand)
        L = sp.expand(M.trace())
        F = 0
        for _ in range(400):
            if L == 0:
                break
            P = sp.Poly(sp.expand(L * t ** 400), t)
            d = P.degree() - 400
            if d < 0:
                raise ValueError('trace is not symmetric in t, 1/t')
            c = P.coeff_monomial(t ** (d + 400))
            F += c * z ** d
            L = sp.expand(L - c * (t + 1 / t) ** d)
        return sp.expand(F)

    assert (fricke('A'), fricke('B'), fricke('AB')) == (x, y, z)
    E = y ** 2 + z ** 2 - 3
    trb = fricke(W_B)
    check('tr b = yz - x, so tr rho(b) = 0 reads x = yz', sp.expand(trb - (y * z - x)) == 0)
    on = lambda F: sp.expand(F.subs(x, y * z))
    check('on x = yz: tr a = -y E', sp.expand(on(fricke(W_A)) + y * E) == 0)
    check('on x = yz: tr g = -z E', sp.expand(on(fricke(W_G)) + z * E) == 0)
    check('on x = yz: tr(g a) = yz E^2', sp.expand(on(fricke(W_G + W_A)) - y * z * E ** 2) == 0)
    rem_w = sp.rem(sp.Poly(on(fricke(W_W)) - 2, z), sp.Poly(E, z)).as_expr()
    check('tr w - 2 lies in the ideal (x - yz, E)', sp.expand(rem_w) == 0)
    w4 = inv(W_BA) * 2 + W_C + W_BA * 2 + W_A
    rem4 = sp.expand(sp.rem(sp.Poly(on(fricke(w4)) - 2, z), sp.Poly(E, z)).as_expr())
    print(f'      (ba)^-2 c (ba)^2 a: tr - 2 reduces to {sp.factor(rem4)}')
    check('for (ba)^-2 c (ba)^2 a the reduction leaves -y^2 - 1 (not trivial on C_0)',
          sp.expand(rem4 - (-y ** 2 - 1)) == 0)


# ---------------------------------------------------------------------------------------------- exact points
def part4():
    print('(4) rho(w) = 1 at three points of C_0 with algebraic coordinates (exact)')
    s3 = sp.sqrt(3)
    for s in (sp.Rational(1, 3), sp.Rational(2, 7), sp.Rational(-5, 4)):
        cps, sps = (1 - s ** 2) / (1 + s ** 2), 2 * s / (1 + s ** 2)
        Y, Z = s3 * cps, s3 * sps                       # a point of the ellipse y^2 + z^2 = 3 ...
        X = Y * Z                                       # ... with x = yz, so E = 0 and tr rho(b) = 0
        cu, cv = X / 2, Y / 2
        su, sv = sp.sqrt(1 - cu ** 2), sp.sqrt(1 - cv ** 2)
        tau = (cu * cv - Z / 2) / (su * sv)             # tr(AB) = 2 (cos u cos v - tau sin u sin v) = Z
        nt = sp.sqrt(1 - tau ** 2)
        word = {'A': (cu, su, 0, 0), 'a': (cu, -su, 0, 0),
                'B': (cv, sv * tau, sv * nt, 0), 'b': (cv, -sv * tau, -sv * nt, 0)}
        r = [1, 0, 0, 0]
        for ch in W_W:
            r = [sp.expand(c) for c in qm(r, word[ch])]
        r = [sp.simplify(c) for c in r]
        inside = abs(sp.N(tau, 30)) < 1 and Y != 0
        check(f's = {s}: (y, z) = ({sp.nsimplify(Y)}, {sp.nsimplify(Z)}), |tau| < 1, y != 0, rho(w) = {r}',
              inside and r == [1, 0, 0, 0])


# ---------------------------------------------------------------------------------------------- 50 digits
def part5():
    print('(5) 50-digit evaluation on the parametrization of C_0 (Lemma 3.3 and Definition 3.4 of the companion paper)')
    mp.mp.dps = 50

    def mq(p, q):
        w1, x1, y1, z1 = p
        w2, x2, y2, z2 = q
        return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2)

    def qe(t, ax):
        return (mp.cos(t), mp.sin(t) * ax[0], mp.sin(t) * ax[1], mp.sin(t) * ax[2])

    def point(v, branch):
        Qv = mp.sin(5 * v) / mp.sin(v)
        u = mp.acos(branch * mp.sqrt(max(mp.mpf(0), 1 - Qv)) / 2)
        tau = mp.cot(2 * v) * mp.cot(u)
        QA, QB = (1, 0, 0), (tau, mp.sqrt(1 - tau ** 2), 0)
        rA = lambda e: qe(e * u, QA)
        rB = lambda e: qe(e * v, QB)
        a = mq(rA(2), rB(3))
        b = mq(rB(-2), rA(-1))
        c = mq(mq(rB(-2), a), rB(2))
        av, bv, cv = a[1:], b[1:], c[1:]
        dot = lambda p, q: sum(pi * qi for pi, qi in zip(p, q))
        cg = dot(av, bv)
        g = mp.acos(cg)
        e2 = [(bi - cg * ai) / mp.sin(g) for ai, bi in zip(av, bv)]
        th = mp.atan2(dot(cv, e2), dot(cv, av))
        return tau, g, th

    worst, count = mp.mpf(0), 0
    for branch in (1, -1):
        for i in range(1, 400):
            v = mp.pi / 6 + (2 * mp.pi / 3) * i / 400
            if abs(v - mp.pi / 2) < mp.mpf('1e-30'):
                continue
            tau, g, th = point(v, branch)
            if abs(tau) >= 1:
                continue
            dev = mp.fmod(th - 6 * g + mp.pi, 2 * mp.pi)
            worst = max(worst, min(abs(dev), abs(abs(dev) - 2 * mp.pi)))
            count += 1
    print(f'      {count} points, max |theta - 6 gamma + pi| mod 2 pi = {mp.nstr(worst, 5)}')
    check('|theta - 6 gamma + pi| <= 1.4e-47 on C_0', worst <= mp.mpf('1.4e-47'))
    for half, (lo, hi) in (('L', (mp.pi / 6, mp.pi / 2)), ('R', (mp.pi / 2, 5 * mp.pi / 6))):
        vs = [lo + (hi - lo) * i / 600 for i in range(1, 600)]
        order = [(1, list(reversed(vs))), (-1, vs)] if half == 'L' else [(1, vs), (-1, list(reversed(vs)))]
        seq = []
        for br, vv in order:
            for v in vv:
                tau, g, th = point(v, br)
                if abs(tau) < 1:
                    seq.append(g)
        d = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
        mono = all(x > 0 for x in d) or all(x < 0 for x in d)
        check(f'half C_0^{half}: gamma runs from {mp.nstr(seq[0] / mp.pi, 6)} pi to {mp.nstr(seq[-1] / mp.pi, 6)} pi, '
              f'strictly monotone', mono)


def main():
    t0 = time.time()
    part1()
    part3()
    part4()
    part5()
    ok = all(checks)
    print(f"{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}: {sum(checks)}/{len(checks)} ({time.time() - t0:.1f} s)")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
