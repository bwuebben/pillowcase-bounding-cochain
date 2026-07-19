#!/usr/bin/env python3
"""
Fintushel-Stern / Poudel-Saveliev Z/4 gradings of the singular-instanton CHAIN
complex IC^natural(T(3,n)), n odd.

Self-computes the individual spectral-flow gradings mu(alpha) of the irreducible
flat connections on the double branched cover Sigma(2,3,n), and assembles the
Z/4-graded singular instanton CHAIN complex IC^natural = (1+a,a,a,a).

NOTE (correction 2026-07-14): this is the CHAIN complex, not the homology. The
homology I^natural equals the chain (differential zero) only for n = 1 (mod 6);
for n = 5 (mod 6) the differential has rank 1 and rank I^natural = 4a-1 = (1+4a)-2
(e.g. T(3,5) = P(-2,3,5) = 10_124 has rank 7, not 9). See RESEARCH_LOG.md sec 19.

Verified against Anvari (arXiv:1609.05025) Example 6.1:
  Sigma(2,3,7): (1,2,2) -> gr=175 (=7 mod 8), mu=87 (=3 mod 4)
                (1,2,4) -> gr=251 (=3 mod 8), mu=125 (=1 mod 4)

Formula (p=3, q=n):
  gr(alpha) = e^2/(pq) + (2/p) S_p + (2/q) S_q,   e = pq*l1 + 2q*l2 + 2p*l3,
    S_r = sum_{k=1}^{r-1} cot(pi k/r) cot(pi b_r k/r) sin^2(pi k l_r/r),
    b_2 = (2q)^{-1} mod p,  b_3 = (2p)^{-1} mod q     (Seifert framing inverses).
  rho(alpha) = 1 - 2/(pq) - (4/p) sum cot·tan·cos^2 - (4/q) sum cot·tan·cos^2  (Anvari Thm B).
  mu(alpha) = (1/2) gr + (1/4)(1 - rho)   (mod 4).
Each irreducible contributes 4 generators (2 at mu, 2 at mu+1); theta -> 1 at sigma mod 4 (=0 for n odd).
"""
import math
from math import pi, cos, sin, tan

def cot(x): return cos(x) / sin(x)

def inv_mod(x, m):
    x %= m
    for t in range(1, m):
        if (x * t) % m == 1:
            return t
    return None  # not invertible (n even => Sigma not a homology sphere; formula N/A)

def gr_FS(q, l1, l2, l3, p=3):
    b2, b3 = inv_mod(2 * q, p), inv_mod(2 * p, q)
    e = p * q * l1 + 2 * q * l2 + 2 * p * l3
    Sp = sum(cot(pi * k / p) * cot(pi * b2 * k / p) * sin(pi * k * l2 / p) ** 2 for k in range(1, p))
    Sq = sum(cot(pi * k / q) * cot(pi * b3 * k / q) * sin(pi * k * l3 / q) ** 2 for k in range(1, q))
    return round(e * e / (p * q) + (2.0 / p) * Sp + (2.0 / q) * Sq)

def rho_FS(q, l2, l3, p=3):
    b2, b3 = inv_mod(2 * q, p), inv_mod(2 * p, q)
    r = 1 - 2.0 / (p * q)
    r -= (4.0 / p) * sum(cot(pi * k / p) * tan(pi * b2 * k / p) * cos(pi * k * l2 / p) ** 2 for k in range(1, p))
    r -= (4.0 / q) * sum(cot(pi * k / q) * tan(pi * b3 * k / q) * cos(pi * k * l3 / q) ** 2 for k in range(1, q))
    return r

def mu(q, l2, l3):
    g, r = gr_FS(q, 1, l2, l3), rho_FS(q, l2, l3)
    return round(0.5 * g + 0.25 * (1 - r)) % 4, g, r

def valid(q, l3, p=3):  # spherical-triangle inequalities for (1,2,l3) on Sigma(2,3,q)
    ang = [pi / 2, 2 * pi / 3, pi * l3 / q]
    return (ang[0] < ang[1] + ang[2] and ang[1] < ang[0] + ang[2]
            and ang[2] < ang[0] + ang[1] and sum(ang) < 2 * pi)

if __name__ == "__main__":
    print("Verify Sigma(2,3,7) vs Anvari Ex 6.1:")
    for l3 in (2, 4):
        m, g, r = mu(7, 2, l3)
        print(f"  (1,2,{l3}): gr={g} rho={r:.3f} mu={m}")
    print("\nT(3,n), n odd: self-computed gradings and assembled CHAIN complex IC^natural = (1+a,a,a,a)")
    print("  (homology I^natural = chain for n=1 mod 6; for n=5 mod 6 rank drops by 2 to 4a-1):")
    for q in range(5, 44, 2):
        if q % 3 == 0:
            continue
        irr = [l3 for l3 in range(2, q, 2) if valid(q, l3)]
        mus = sorted(mu(q, 2, l3)[0] for l3 in irr)
        C = [0, 0, 0, 0]
        for m in mus:
            C[m] += 2
            C[(m + 1) % 4] += 2
        C[0] += 1  # theta, sigma = -4a = 0 mod 4
        a = len(irr)
        assert C == [1 + a, a, a, a] and (C[0] - C[1] + C[2] - C[3]) == 1
        homology = (1 + 4 * a) if q % 6 == 1 else (4 * a - 1)  # differential fires iff n=5 mod 6
        print(f"  T(3,{q:2d}) [n%6={q%6}]: a={a:2d}  mu={mus}  IC^nat={tuple(C)} (rank {1+4*a})"
              f"  ->  rank I^nat={homology}")
