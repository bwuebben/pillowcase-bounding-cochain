# Programs for Sections 4 and 5 and Appendix A

These programs carry out the computations of Sections 4 and 5 and Appendix A of
*Instanton and pillowcase homology of the (−2,3,q) pretzel knots*. Statement
numbers refer to the compiled paper.

## The sheared decompositions

These scripts reproduce the computations reported in Sections 4 and 5 and in Appendix A.1(i)–(v). Statement numbers
refer to the compiled paper. The proofs of Theorems 1.2 and 1.3 (Theorems 4.18 and 5.1) are by hand, except for the base case of
Theorem 1.2 (Remark 4.3), which is the verified run of `pretzel_chirality_snappy.py`; none of the other
computations is used in them. There are two kinds of computation.

* **Exact computations.** Symbolic identities (sympy), exact rational arithmetic and exact free-group words. They check,
  step by step or over a stated range, the displayed identities and the combinatorial steps of the hand proofs: the
  straightness of Lemma 4.7, the crossing combinatorics of Propositions 4.13 and 4.17, the group identities of
  Lemmas 4.2 and 5.7, and the Alexander polynomials, determinants and signatures of Propositions 5.4 and 5.5. One
  further step is a 50-digit evaluation.
* **Floating-point numerics.** The perturbed character varieties are traced by numerical continuation, and the
  complexes of the sheared curves are computed from the traced curves (Remarks 4.21 and 5.8, Appendix A.1(iv)–(v)).
  They corroborate Theorems 4.18 and 5.1 in examples and are not part of the proofs.

Each script prints the numbers of the statements it checks next to the computed values, with `PASS` or `FAIL`, and exits with
status 0 if and only if all its checks pass.

| script | verifies | runtime |
|---|---|---|
| `straightness_identity.py` | Lemma 4.7 and Appendix A.1(i), exact: every displayed identity of Steps 1–7 of the proof, symbolically in (v, ζ). Without E = 0: vec ρ(g)·vec ρ(a) = ½ζ sin 2v E² and tr(ρ(g)ρ(a)) = −2ζ sin 2v E². Fricke form: tr b = yz − x; on x = yz, tr a = −yE, tr g = −zE, tr(ga) = yzE², tr w − 2 ∈ (x − yz, E); for (ba)⁻²c(ba)²a the remainder is −y² − 1. ρ(w) = 1 at three algebraic points of C₀. At 50 digits, \|θ − 6γ + π\| ≤ 1.4·10⁻⁴⁷, and γ is strictly monotone on each half of C₀ | 7 s |
| `sheared_polygon_counts.py` | Appendix A.1(ii), exact, for every odd 7 ≤ q ≤ 401. The values (4.5) and the margins of Remark 4.12. The crossing lists of Proposition 4.13 for both resolutions: N_q triple lines, q + 2 + 4N_q generators on the arc (joined case), q − 2 on the arc and 4 + 4N_q on the circle (split case), (q + 2)π − φ₂ = (q + 18)π/6, no crossing at γ = π/2. The abscissae, differences (≥ π/(q(q − 6))) and order of Proposition 4.17. The degrees of Proposition 4.15 and the graded homology (1 + N₃, N₁, N₁, N₃) with rank 2N_q of Theorem 4.18(3). The identity of Remark 4.19. Also the values N₇, …, N₁₉ of Section 4.5 | 1 s |
| `pretzel_groups.py` | Exact free-group words. ba = cd for the words of the tangles of Sections 4–5 (Appendix A.1(iii)). Lemma 4.2 for n ∈ {5, 8, …, 29}, m ≤ 15: the relator R_{n,m}, the Clay–Watson relator (4.4) under α = A, β = B⁻¹, and the control that t_e^{+m} never gives (4.4). Lemma 5.7(b) as a word identity, ν, m ≤ 15. Fox calculus against Burau for all 4 ≤ n ≤ 26 prime to 3 and m ≤ 6 (112 pairs) (Appendix A.1(iii)) | 5 s |
| `twisted_alexander.py` | Propositions 5.4 and 5.5 and Lemmas 5.3 and 5.6, exact. The Burau matrix facts and the telescoping identities of the proof of Proposition 5.4, as identities of rational functions, and formula (5.1). The closed form for all n ≤ 80 prime to 3 and m ≤ 20 (1092 knots), with coefficients in {−1, 0, 1}, degree 2(n + m − 1), ℓ = 4ν + 2M + 1 = \|σ\| + 1. The example T(3,7;2,2). The matrices at t = −1 in the proof of Proposition 5.5, and det = 2m + 1 or 2m + 3 for 7 ≤ n ≤ 50, m ≤ 20. The case n = 5 against (3.2). The inputs of Lemma 5.6 | 7 s |
| `sheared_shear_check.py` | Lemma 4.6(a) and Appendix A.1(iv), numerical: on the traced variety, t_e^h acts by S_h to within 10⁻¹² for h = −6, …, 2, while the opposite shear is off by about π | 1 s |
| `sheared_resolution_type.py` | Proposition 4.10 and Appendix A.1(v): dκ_c(0) = (−τ, −½) at c±, from the derivative routine of the continuation code. The type of resolution at \|ε\| = 0.03 in 24 directions: four excluded, agreement in the other 20, split at 75°, 90°, 105°, 255°, 270°, 285° | 18 s |
| `sheared_numerics.py` | Appendix A.1(v) (summarized in Remark 4.21), joined sectors, numerical: ten perturbations, q ∈ {7, 9, …, 25, 29, 33}, earrings 0.002 and 0.0005 (240 cases). Graded homology as in Theorem 4.18 in all 240 cases. Bigons from down- to up-crossings with #bigons = 2 rank d. Degrees as in Proposition 4.15. The generator count q + 2 + 4N_q in 178 cases and short in 62, only at the larger perturbations. `--quick`: 12 cases | 136 s (`--quick` 10 s) |
| `sheared_split.py` | Appendix A.1(v) (summarized in Remark 4.21), split sectors, numerical: five perturbations, q = 7, …, 21 (40 cases). The arc is φ-monotone with q − 2 generators. The homology of the remnant circle, computed from its own complex, has rank 1 in each degree. The total is as in Theorem 4.18. The circle carries 4 + 4N_q generators in 37 cases, the exceptions being q = 17 at three perturbations. A lift of L₀ meets the lift of the circle up to seven times (q = 21) | 19 s |
| `twisted_numerics.py` | Remark 5.8 and Appendix A.1(v), numerical: T(3,n;2,m) for n = 7, 8, 13, m = 1, 2, 3, at two perturbations (18 cases). The arc carries 1 + 2M generators and no bigon. There are 2j circles, each with homology (1,1,1,1) computed from its complex. The total is the vector of Theorem 5.1(iii) | 34 s |
| `pretzel_chirality_snappy.py` | Remark 4.3 (needs SnapPy; `--verified` needs SnapPy inside Sage and lists all isometries with interval arithmetic: two, both preserving orientation and meridians, with a Rolfsen-twist check of the surgery convention). From a planar diagram of [HHK2, Fig. 18], the exterior of K⁽⁰⁾ ∪ e is isometric to that of the closure of the positive braid (σ₁σ₂)⁵ with a circle around two strands. The isometries preserve orientation and meridians and extend to the links, and an exact combinatorial isomorphism of canonical triangulations is found. Also, the fillings of e along (1, h), h = −8, …, 1, have the Alexander polynomials of P(−2,3,5 − 2h). With `--verified` this is the proof of Remark 4.3 | < 1 s |

Shared code:

* `t3n_pillowcase.py`: the continuation code of the companion paper on T(3,n), copied unchanged from its code
  directory. It provides the equations of the tangle of Hedden, Herald and Kirk, the restriction map, the tracing of
  W_ε, and that paper's complex code (used here only by `twisted_numerics.py`).
* `sheared_complex.py`: independent code for the sheared curves of Section 4. It covers the restriction map from the
  quaternion representation and the lift, the shear S_{−m}, and the lifts of L₀ and the generators. Degrees use
  Maslov indices computed on the unsheared curve against the line field spanned by (1, 6 − q) (Lemma 4.6(f)); an index
  obtained by unwrapping tangent angles on the sheared polygonal curve is unreliable near the junctions. It also
  enumerates bigons as embedded discs in ℝ² \ (πℤ)² and computes graded homology over 𝔽₂.

## Requirements and use

Python 3.11 or later, with numpy, sympy and mpmath (`requirements.txt`). SnapPy is needed only for
`pretzel_chirality_snappy.py`, which is reported as `SKIP` without it. Tested with Python 3.11.14 and 3.14.4
(numpy 2.4.6 and 2.5.3, sympy 1.14.0, mpmath 1.3.0, snappy 3.3.2, spherogram 2.4.1).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PY=.venv/bin/python ./run_all.sh            # everything, about 4 minutes
PY=.venv/bin/python ./run_all.sh --quick    # shorter run of sheared_numerics.py, about 2 minutes
```
