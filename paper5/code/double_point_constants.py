#!/usr/bin/env python3
r"""Numerical values stated after Lemma 4.13 and after Lemma 5.4, at the double points c_+, c_- of the
unperturbed variety W (n = 4, 5 mod 6).

 1. After Lemma 4.13: the leading coefficients V_gamma, V_phi of (gamma, phi) - (gamma(c), phi(c)) = +-(V_gamma, V_phi) w^2
    along the central circle, w = |v - pi/2|.  For n = 4: V_gamma = 2/sqrt 3, V_phi = 2 sqrt 3, slope 4.  For
    n = 4, 5, 10, 11, 16, 17, 41 the coefficient of phi given by the formula is compared with the one obtained by
    evaluating the restriction map along the central circle (agreement to at least four significant digits).
 2. After Lemma 5.4: the first-order coefficients kappa_c = lambda_c eps_A + mu_c eps_B + O(|eps|^2) of the critical
    value at the double points, with the covector ell_c normalized to unit length: mu_c = +-1/4 for n = 4, 5;
    |mu_c| = 0.1179, 0.0772, 0.0574 for n = 10, 16, 22 and for n = 11, 17, 23; lambda_c = 0 for n = 4 (mod 6)
    and |lambda_c| ~ 0.44 for n = 5 (mod 6).

Floating-point numerics.
"""
import math
import sys

import numpy as np

from t3n_pillowcase import PI, TAU, Tangle, decomposition

results = []


def check(label, ok, computed, paper):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: computed {computed}; paper {paper}")


def params(n):
    eta = 1 if n % 3 == 1 else -1
    m = (n - eta) // 3
    return eta, m, 2 * m + eta, m + eta                    # eta, m, m', k


def gamma_star(n):
    eta, m, mp, k = params(n)
    return math.acos(math.sqrt((2 * m + 1) / (2 * m + 2))) if n % 6 == 4 else math.acos(math.sqrt(1 - 1 / (2 * m)))


def V_formula(n):
    """(V_gamma, V_phi) of Lemma 4.13."""
    eta, m, mp, k = params(n)
    h = 1 / m
    if n % 6 == 4:
        kappa = (h + 2) * (3 * h + 1) / (6 * (h + 1))
        c3 = 4 * (h + 1) * (h + 2) * (h + 3)
        Vphi = c3 * m ** 3 / (8 * math.sqrt(mp) * k)
    else:
        kappa = (2 - h) * (2 * h + 1) / 6
        c3 = 4 / 3 * (2 - h) * (11 - 2 * h)
        Vphi = c3 * m ** 2 / (8 * math.sqrt(mp))
    return kappa * h * m ** 2 / math.sin(2 * gamma_star(n)), Vphi


def central_point(n, v, branch):
    """the point of the central circle over v on the branch sign(cos u) = branch (Lemma 3.3)."""
    eta, m, mp, k = params(n)
    Q = math.sin(n * v) / math.sin(k * v)
    x = branch * math.sqrt((1 - Q) / 4)
    u = math.acos(x)
    return np.array([u, v, math.cos(m * v) / math.sin(m * v) * x / math.sin(u)])


def wrap(a):
    return (a + PI) % TAU - PI


def V_direct(n):
    """leading coefficients in w^2 of gamma and phi along both branches of the left half near v = pi/2,
    from a polynomial fit of the restriction map; returns {('+' or '-'): (Vg, Vphi)} keyed by the crossing."""
    r, s = decomposition(n)
    T = Tangle(3, n, r, s)
    gs = gamma_star(n)
    out = {}
    ws = np.geomspace(1e-3, 1e-2, 40) / max(1, params(n)[1] / 4)
    for branch in (1, -1):
        P = np.array([T.pillow(central_point(n, PI / 2 - w, branch)) for w in ws])
        g0 = P[0, 0]
        cross = '+' if abs(g0 - gs) < abs(g0 - (PI - gs)) else '-'
        gc = gs if cross == '+' else PI - gs
        phic = (-3 * gc) if n % 6 == 4 else -gc               # phi on the image of B (Proposition 3.9)
        dg = P[:, 0] - gc
        dphi = wrap((P[:, 1] - P[:, 0]) - phic)
        Vand = np.vander(ws, 8, increasing=True)[:, 2:]          # w^2 ... w^7
        a = np.linalg.lstsq(Vand, dg, rcond=None)[0][0]
        b = np.linalg.lstsq(Vand, dphi, rcond=None)[0][0]
        out[cross] = (a, b)
    return out


def main():
    print("After Lemma 4.13: leading coefficients at the double points (floating-point numerics)")
    Vg, Vp = V_formula(4)
    check("n = 4: V_gamma = 2/sqrt(3), V_phi = 2 sqrt(3), slope 1 + V_phi/V_gamma = 4",
          abs(Vg - 2 / math.sqrt(3)) < 1e-12 and abs(Vp - 2 * math.sqrt(3)) < 1e-12 and abs(1 + Vp / Vg - 4) < 1e-12,
          f"V_gamma = {Vg:.6f}, V_phi = {Vp:.6f}, slope {1 + Vp / Vg:.6f}", "2/sqrt3, 2 sqrt3, 4")
    for n in (4, 5, 10, 11, 16, 17, 41):
        Vg, Vp = V_formula(n)
        D = V_direct(n)
        for cross, sgn in (('+', 1), ('-', -1)):
            dg, dp = D[cross]
            rel = abs(dp - sgn * Vp) / Vp
            relg = abs(dg - sgn * Vg) / Vg
            check(f"n = {n:2d}, c_{cross}: leading coefficient of phi (sign {'+' if sgn > 0 else '-'})", rel < 5e-4,
                  f"direct {dp:+.6f}, formula {sgn * Vp:+.6f} (rel. diff. {rel:.1e}); gamma: direct {dg:+.6f}, "
                  f"formula {sgn * Vg:+.6f}", "agreement to >= 4 significant digits")
    print("After Lemma 5.4: first-order coefficients of kappa_c (unit covector ell_c)")
    paper_mu = {4: 0.25, 5: 0.25, 10: 0.1179, 16: 0.0772, 22: 0.0574, 11: 0.1179, 17: 0.0772, 23: 0.0574}
    for n in (4, 10, 16, 22, 5, 11, 17, 23):
        eta, m, mp, k = params(n)
        r, s = decomposition(n)
        T = Tangle(3, n, r, s)
        if n % 6 == 4:
            u1 = math.acos(math.sqrt((2 * m + 1) / (2 * m + 2)))
            cs = [np.array([u1, PI / 2, 0.0]), np.array([PI - u1, PI / 2, 0.0])]
        else:
            t1 = math.sqrt(1 - 1 / (2 * m))
            cs = [np.array([PI / 2, PI / 2, t1]), np.array([PI / 2, PI / 2, -t1])]
        vals = []
        for c in cs:
            U, sv, _ = np.linalg.svd(T.jac(c))
            ell = U[:, -1]
            lam, mu = ell @ T.dpsi_deps(c)
            vals.append((sv[-1], lam, mu))
        rank1 = all(v[0] < 1e-12 and abs(T.psi(c)).max() < 1e-12 for v, c in zip(vals, cs))
        mus = [abs(v[2]) for v in vals]
        lams = [abs(v[1]) for v in vals]
        ok = rank1 and all(round(x, 4) == paper_mu[n] for x in mus)
        if n % 6 == 4:
            ok &= all(x < 1e-12 for x in lams)
            lam_txt, lam_paper = "lambda_c = 0", "lambda_c = 0"
        else:
            ok &= all(abs(x - 0.44) < 0.01 for x in lams)
            lam_txt, lam_paper = f"|lambda_c| = {lams[0]:.4f}, {lams[1]:.4f}", "|lambda_c| ~ 0.44"
        check(f"n = {n:2d}: d Psi has rank 1 at c_+-; |mu_c|, lambda_c", ok,
              f"|mu_c| = {mus[0]:.6f}, {mus[1]:.6f}; {lam_txt}", f"|mu_c| = {paper_mu[n]}; {lam_paper}")
    # Remark 5.8: the sampled perturbations satisfy |lambda_c eps_A| > 2 |mu_c eps_B| at both double points for n = 5 (mod 6)
    samples = ((0.5, 0.03), (0.3, 0.05), (-0.5, 0.03), (0.5, -0.03), (-0.3, 0.03))   # (n eps_A, eps_B)
    ratios = []
    for n in range(5, 48, 6):
        eta, m, mp, k = params(n)
        r, s = decomposition(n)
        T = Tangle(3, n, r, s)
        t1 = math.sqrt(1 - 1 / (2 * m))
        for c in (np.array([PI / 2, PI / 2, t1]), np.array([PI / 2, PI / 2, -t1])):
            U, sv, _ = np.linalg.svd(T.jac(c))
            lam, mu = U[:, -1] @ T.dpsi_deps(c)
            ratios += [(abs(lam * nA / n) / abs(mu * eB), n) for nA, eB in samples]
    worst = min(ratios)
    check("n = 5, 11, ..., 47 at the sampled perturbations: |lambda_c eps_A| > 2 |mu_c eps_B|", worst[0] > 2,
          f"minimum ratio {worst[0]:.4f} (n = {worst[1]})", "> 2")
    print("\nConstants at the double points: %d checks, %d passed" % (len(results), sum(results)))
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
