# Programs for Sections 4, 5 and 7 and Appendix A

These programs carry out the computations of Sections 4, 5 and 7, Remark 6.17 and Appendix A of
*Instanton and pillowcase homology of the (−2,3,q) pretzel knots*. Statement numbers refer to the
compiled paper. The programs for the rest of Section 6 are in the directory `pillowcase/`.

Every script prints each number the paper states next to the computed value, with `PASS` or
`FAIL`, and exits with status 0 if and only if all its checks pass.

## Sections 4 and 5: the sheared decompositions

These scripts reproduce the computations reported in Sections 4 and 5 and in Appendix A.3(i)–(v). Statement numbers
refer to the compiled paper. The proofs of Theorems 1.2 and 1.3 (Theorems 4.18 and 5.1) are by hand; none of these
computations is used in them. There are two kinds of computation.

* **Exact computations.** Symbolic identities (sympy), exact rational arithmetic and exact free-group words. They check,
  step by step or over a stated range, the displayed identities and the combinatorial steps of the hand proofs: the
  straightness of Lemma 4.7, the crossing combinatorics of Propositions 4.13 and 4.17, the group identities of
  Lemmas 4.2 and 5.7, and the Alexander polynomials, determinants and signatures of Propositions 5.4 and 5.5. One
  further step is a 50-digit evaluation.
* **Floating-point numerics.** The perturbed character varieties are traced by numerical continuation, and the
  complexes of the sheared curves are computed from the traced curves (Remarks 4.21 and 5.8, Appendix A.3(iv)–(v)).
  They corroborate Theorems 4.18 and 5.1 in examples and are not part of the proofs.

Every script prints each number the paper states next to the computed value, with `PASS` or `FAIL`, and exits with
status 0 if and only if all its checks pass.

| script | verifies | runtime |
|---|---|---|
| `straightness_identity.py` | Lemma 4.7 and Appendix A.3(i), exact: every displayed identity of Steps 1–7 of the proof, symbolically in (v, ζ). Without E = 0: vec ρ(g)·vec ρ(a) = ½ζ sin 2v E² and tr(ρ(g)ρ(a)) = −2ζ sin 2v E². Fricke form: tr b = yz − x; on x = yz, tr a = −yE, tr g = −zE, tr(ga) = yzE², tr w − 2 ∈ (x − yz, E); for (ba)⁻²c(ba)²a the remainder is −y² − 1. ρ(w) = 1 at three algebraic points of C₀. At 50 digits, \|θ − 6γ + π\| ≤ 1.4·10⁻⁴⁷, and γ is strictly monotone on each half of C₀ | 7 s |
| `sheared_polygon_counts.py` | Appendix A.3(ii), exact, for every odd 7 ≤ q ≤ 401. The values (4.5) and the margins of Remark 4.12. The crossing lists of Proposition 4.13 for both resolutions: N_q triple lines, q + 2 + 4N_q generators on the arc (joined case), q − 2 on the arc and 4 + 4N_q on the circle (split case), (q + 2)π − φ₂ = (q + 18)π/6, no crossing at γ = π/2. The abscissae, differences (≥ π/(q(q − 6))) and order of Proposition 4.17. The degrees of Proposition 4.15 and the graded homology (1 + N₃, N₁, N₁, N₃) with rank 2N_q of Theorem 4.18(3). The identity of Remark 4.19. Also the values N₇, …, N₁₉ of Section 4.5 | 1 s |
| `pretzel_groups.py` | Exact free-group words. ba = cd for the words of the tangles of Sections 4–5 (Appendix A.3(iii)). Lemma 4.2 for n ∈ {5, 8, …, 29}, m ≤ 15: the relator R_{n,m}, the Clay–Watson relator (4.4) under α = A, β = B⁻¹, and the control that t_e^{+m} never gives (4.4). Lemma 5.7(b) as a word identity, ν, m ≤ 15. Fox calculus against Burau for all 4 ≤ n ≤ 26 prime to 3 and m ≤ 6 (112 pairs) (Appendix A.3(iii)) | 5 s |
| `twisted_alexander.py` | Propositions 5.4 and 5.5 and Lemmas 5.3 and 5.6, exact. The Burau matrix facts and the telescoping identities of the proof of Proposition 5.4, as identities of rational functions, and formula (5.1). The closed form for all n ≤ 80 prime to 3 and m ≤ 20 (1092 knots), with coefficients in {−1, 0, 1}, degree 2(n + m − 1), ℓ = 4ν + 2M + 1 = \|σ\| + 1. The example T(3,7;2,2). The matrices at t = −1 in the proof of Proposition 5.5, and det = 2m + 1 or 2m + 3 for 7 ≤ n ≤ 50, m ≤ 20. The case n = 5 against (3.2). The inputs of Lemma 5.6 | 7 s |
| `sheared_shear_check.py` | Lemma 4.6(a) and Appendix A.3(iv), numerical: on the traced variety, t_e^h acts by S_h to within 10⁻¹² for h = −6, …, 2, while the opposite shear is off by about π | 1 s |
| `sheared_resolution_type.py` | Proposition 4.10 and Remark 4.21: dκ_c(0) = (−τ, −½) at c±, from the derivative routine of the continuation code. The type of resolution at \|ε\| = 0.03 in 24 directions: four excluded, agreement in the other 20, split at 75°, 90°, 105°, 255°, 270°, 285° | 18 s |
| `sheared_numerics.py` | Remark 4.21, joined sectors, numerical: ten perturbations, q ∈ {7, 9, …, 25, 29, 33}, earrings 0.002 and 0.0005 (240 cases). Graded homology as in Theorem 4.18 in all 240 cases. Bigons from down- to up-crossings with #bigons = 2 rank d. Degrees as in Proposition 4.15. The generator count q + 2 + 4N_q in 178 cases and short in 62, only at the larger perturbations. `--quick`: 12 cases | 136 s (`--quick` 10 s) |
| `sheared_split.py` | Remark 4.21, split sectors, numerical: five perturbations, q = 7, …, 21 (40 cases). The arc is φ-monotone with q − 2 generators. The homology of the remnant circle, computed from its own complex, has rank 1 in each degree. The total is as in Theorem 4.18. The circle carries 4 + 4N_q generators in 37 cases, the exceptions being q = 17 at three perturbations. A strand meets the lift of the circle up to seven times (q = 21) | 19 s |
| `twisted_numerics.py` | Remark 5.8 and Appendix A.3(v), numerical: T(3,n;2,m) for n = 7, 8, 13, m = 1, 2, 3, at two perturbations (18 cases). The arc carries 1 + 2M generators and no bigon. There are 2j circles, each with homology (1,1,1,1) computed from its complex. The total is the vector of Theorem 5.1(iii) | 34 s |
| `pretzel_chirality_snappy.py` | Remark 4.3 (needs SnapPy). From a planar diagram of [HHK2, Fig. 18], the exterior of K⁽⁰⁾ ∪ e is isometric to that of the closure of the positive braid (σ₁σ₂)⁵ with a circle around two strands. The isometries preserve orientation and meridians and extend to the links, and an exact combinatorial isomorphism of canonical triangulations is found. Also, the fillings of e along (1, h), h = −8, …, 1, have the Alexander polynomials of P(−2,3,5 − 2h). Not used in the proofs | < 1 s |

Shared code:

* `t3n_pillowcase.py`: the continuation code of the companion paper on T(3,n), copied unchanged from its code
  directory. It provides the equations of the tangle of Hedden, Herald and Kirk, the restriction map, the tracing of
  W_ε, and that paper's complex code (used here only by `twisted_numerics.py`).
* `sheared_complex.py`: independent code for the sheared curves of Section 4. It covers the restriction map from the
  quaternion representation and the lift, the shear S_{−m}, and the strands of L₀ and the generators. Degrees use
  Maslov indices computed on the unsheared curve against the line field spanned by (1, 6 − q) (Lemma 4.6(f)); an index
  obtained by unwrapping tangent angles on the sheared polygonal curve is unreliable near the junctions. It also
  enumerates bigons as embedded discs in ℝ² \ (πℤ)² and computes graded homology over 𝔽₂.

## Section 7, Remark 6.17 and Appendix A: Smith's decomposition and the Khovanov arc

These scripts reproduce the computations reported in Section 7, in Remark 6.17 and in Appendix A.1–A.3, and recheck
the morphism dimensions of Propositions 6.7 and 6.8 with two independent implementations. Statement numbers refer to
the compiled paper. There are two kinds of computation.

* **Exact computations over F₂.** Twisted complexes over the two-arc algebra of Kotelskiy, Watson and Zibrowius are
  finite matrices, the Maurer–Cartan equation is the finite identity (δ + b)² = 0, strict isomorphism of loop-type
  complexes is decided by comparing the words of their components, and morphism homology is computed in full, without
  truncating word length. These computations are the proofs of Theorem 7.1 and of the table in the proof of
  Corollary 7.6, and they establish the statements of Remarks 7.2, 7.3, 7.7, 7.8, 6.17, Proposition 7.10 and
  Appendix A.1–A.2. Their inputs are the complexes printed in the paper, curves encoded by `kwz_encode.py`, the
  invariants of Kotelskiy, Watson and Zibrowius computed by Zibrowius's program kht++, and ranks of reduced Khovanov
  homology computed with Khoca (both shipped as data, see below). The searches of Remarks 7.2 and 7.11 are exhaustive
  within the classes they state. Theorem 7.5 is proved in the paper by a geometric argument; the census of linear
  curves in `kwz_linear_pairings.py` supports Remark 7.7.
* **Floating-point input.** The encodings of the numerical traces of Smith's perturbed curve (Appendix A.3, "Traced
  curves") come from floating-point computations, and the statements about them (Remark 7.11, first sentences, and
  the cases q = 17, 19 of its searches) are corroboration only.

Every script prints each number the paper states next to the computed value, with `PASS` or `FAIL` (checks of
consistency that the paper does not state are marked `expected` instead of `paper`), and exits with status 0 if and
only if all its checks pass.

| script | verifies | runtime |
|---|---|---|
| `kwz_encoder_check.py` | Appendix A.3, "Encodings": for each of the 216 rational tangles Q_s computed by kht++, BN~(Q_s) is the straight lattice arc and Kh~(Q_s) the figure-eight around it. Section 7.1: Φ(Kh~(Q_{−r})) is the earring E_r at 49 slopes (and at all 342 slopes of the set of Corollary 7.6). The encodings of E_{−1/2}, E_{−3/4} are the complexes (6.2), (6.3) | 2 s |
| `kwz_hom_check.py` | Propositions 6.7 and 6.8 and Appendix A.3, "Morphism homology": the dimensions 7, 25; 9, 31; 9, 31; 9, 23; 9, 25 by both methods; the data (e, τ, r) = (228, 251, 227), (228, 251, 224), (228, 251, 228) printed after (6.3); f_L = f_M = f_R = 0 and p_{X,i} = U_X; Remark 6.6 (End(N) has f_M = f_R = 1, End(ι_∘, 0) has f_L = f_M = 1); agreement of the two methods on 34 pairings of Section 7 | < 1 s |
| `bn_structure.py` | Theorem 7.1 (a)–(c): (N_q, b) ≅ α_q ⊕ (c, J₂) with components of 23 + 8 and 15 + 8 generators; BN_q is the slide of α_q at (0,0) through two periods, the only one of the 24 slides, for q = 5, 7, …, 21; the Maurer–Cartan identities for κ_q and (N_q, b + κ_q) ≅ BN_q. Remark 7.3: the slide at (π, 0), (N_q, b + κ′_q) ≅ BN~(T_q), and the pairings 13 against 11 (q = 7) and 9 against 7 (q = 5) with E_{−1}. Appendix A.1: N_5, its ends 0 and 22, its pairings 5 and 15, b_5, the components of (N_5, b_5), the arc of (N_5, b_5 + κ_5) and its pairings 7 and 21. Also N_7 and N_5 against the encodings of the piecewise-linear models. Section 7.2: c and (c, J₂) admit a δ-grading but no bigrading, α_q and BN_q are arcs | < 1 s |
| `corner_terms.py` | Remark 7.2 and Appendix A.2: the one-arrow replacements of (N_q, b) (566 candidates for each b at q = 7, 402 at q = 5), exactly two giving BN_q in each case, as tabulated, with the added arrow labelled R at the end over (0,0); κ₇ common to b_{S18} and b_{S25}; the terms for t > 0. Edit distance 2 from N_q to BN_q with two minimal edits, at q = 7 both R-arrows touching 36. The 110 and 58 algebraic switches and the 221 925 and 32 567 sums of at most three, none ≅ BN_q. The census switches at q = 7 against the algebraic switches | 16 s |
| `closures_rational.py` | (7.1) at the 340 (q = 5) and 342 (q = 7) slopes of Corollary 7.6; the three numbers of Corollary 7.6 agree at every slope; its table at E_{−1} and E_{r_q} (7, 30; 11, 46) and the undeformed values 7, 28 and 11, 44; 108 slopes with odd numerator and denominator; N_q below rk Kh~ at 40 and 43 slopes and never above; the values q + 2 at r = −1/2 and 31 at r = −3/4. `--khoca` recomputes the Khovanov ranks with Khoca | 19 s |
| `kwz_linear_pairings.py` | Remark 7.7: 1887 linear curves for each q = 5, 7, …, 13, of which exactly 22 pair differently with BN_q and with α_q ⊕ (c, J₂), all parallel to c (the 19 side patterns on the line through (0,0) and (π,π), the two constant patterns on the line through (π,0) and (0,π), the simple closed curve), none parallel to α_q; the values for c and for the figure-eight around the arc from (0,0) to (π,π). Remark 7.8: 11/13, 22/26, 4, 52/44 (q = 7) and 7/9, 14/18, 36/28 (q = 5). Proposition 7.10. `--quick`: Remark 7.7 for q = 5, 7 only | 2.5 min (`--quick` 49 s) |
| `sign_positive_t.py` | Remark 6.17: for t > 0 the undeformed pairings 7, 25; (9, 31) at S₁₈, S₂₅, (5, 25) at S₆₉ and (7, 25) at S₇₄; the seven switches reaching nine; among the 82 census smoothings, the 82 opposite smoothings and the 16 span elements, (9, 31) only at S₁₈ and S₂₅; the third closure 69 against 103; in the span only b_{S18}, b_{S25} reach nine. The same quantities for t < 0 | 1 s |
| `kwz_traced_curves.py` | Appendix A.3, "Traced curves", and Remark 7.11, first sentences (floating-point input): for q ≤ 13 the traces at a given q have strictly isomorphic encodings, isomorphic to N_q, with 47 and 55 generators at q = 11, 13; the pairing with E_{−1/2} is q (q = 5, 7), q + 4 (q = 11, 13), q + 8 (q = 17, 19); one switch raises it to q + 2 at q = 5, 7 | < 1 s |
| `switches_search.py` | Remark 7.11: edit distances 2, 2, 4, 4 and at least 6 (q = 17, 19, edit distances ≤ 5 excluded); chains of algebraic switches followed by one replacement: one switch at q = 5, 7, none with at most three at q = 11, 13, the first with four at q = 11, ending at an arc plus a closed curve other than α₁₁ ⊕ (c, J₂); switch distance from N_q to α_q ⊕ (c, J₂) 1 at q = 5, 7 and more than 5 at q = 11; the 35 004 complexes (6 331 up to strict isomorphism) obtained from N₁₁ by at most two switches: none has α₁₁ as a summand and none is compatible with the instanton data at the first closure and on the panel of 25 slopes. `--quick`: without the four-switch search, the switch distance at q = 11, the cases q = 17, 19, and with 2000 of the complexes | 30 min (`--quick` 40 s) |

Shared code:

* `kwz_algebra.py`: the algebra, twisted complexes in the notation of the paper, the parser for the output of kht++,
  words and strict isomorphism of loop-type complexes, and the two implementations of the morphism homology: the
  mapping-cone reduction of Lemma 6.5 (the free ranks f_X and the torsion are read off from quotients by powers of
  U_X, and the rank of β_* is computed in such a quotient) and the structure of Mor(X, Y) as a free complex over
  F₂[H], H = D + S² central.
* `kwz_encode.py`: the encoding of curves in the pillowcase as twisted complexes (planar cover, arc system, Smith's
  coordinates with the corner (0, π) at infinity, the translation Φ), straight arcs, figure-eights, linear curves, the
  earrings E_r, KWZ's rational invariants, α_q, c, the doubled word of c, BN_q and the 24 slides.
* `kwz_objects.py`: the complexes printed in the paper (N_7 = N, E, E_*, the switches of Section 6.3, N_5, b_5, κ_q,
  κ′_q), the slope set of Corollary 7.6, and loaders for the data files.
* `switches_common.py`: algebraic switches, chains up to strict isomorphism, and the exact edit distance (tilings of
  the word of the target by disjoint intervals of the component words).

Data (directory `data/`, text files):

* `kht/T<q>yy.kht`, `kht/T<q>yy/cxBNr-c2`, `kht/T<q>yy/cxKhr-c2` (q = 3, 5, …, 21): the input files for the tangles
  T_q = Q_{1/3} + Q_{1/q} and the output files of kht++ with BN~(T_q) and Kh~(T_q).
* `kht_rational.txt`: the input words and the output of kht++ for the 216 rational tangles of `kwz_encoder_check.py`.
* `khovanov_ranks.txt`: rk Kh~(K_r(q); F₂) for the slopes used (q = 5, 7, 9, 11, 13 and the first closures at q = 17,
  19), computed with Khoca 1.5 from planar diagrams built with SnapPy's spherogram.
* `pl_models.txt`: the encodings of the piecewise-linear models L_{q,PL} (q = 5, 7, 11, 13) and of four census
  smoothings, produced by the programs of the directory `pillowcase/`.
* `q7_census.txt` and `q7_smoothings.txt`: the census of the 82 self-intersection orbits of L_{7,PL} with their
  switches, and the encodings of the census and the opposite smoothing at each orbit, produced by the programs of
  `pillowcase/`.
* `traced_curves.txt`: the encodings of the numerical traces of Smith's curve at q = 5, 7, 11, 13, 17, 19 (several
  perturbation parameters and step sizes), traced with `pillowcase/c3_perturbed.py`.

To regenerate the kht++ data, build kht++ from <https://github.com/cbz20/khtpp> and run it on the files
`data/kht/T<q>yy.kht` (and on the words in `data/kht_rational.txt`); the files `cxBNr-c2` and `cxKhr-c2` it writes are
the ones shipped here. To regenerate the Khovanov ranks, install Khoca 1.5 and SnapPy (`pip install khoca snappy`) and
run `python3 closures_rational.py --khoca` (`--sample N` for a random subset).

## Requirements and use

Python 3.11 or later. The programs for Sections 4 and 5 need numpy, sympy and mpmath
(`requirements.txt`); those for Section 7 use only the standard library. SnapPy is needed only for
`pretzel_chirality_snappy.py` (reported as `SKIP` without it), and Khoca and SnapPy only for
`closures_rational.py --khoca`. Tested with Python 3.11.14 and 3.14.4 (numpy 2.4.6 and 2.5.3, sympy 1.14.0,
mpmath 1.3.0, snappy 3.3.2, spherogram 2.4.1).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PY=.venv/bin/python ./run_all.sh            # everything, about 35 minutes
PY=.venv/bin/python ./run_all.sh --quick    # shorter runs of the three longest scripts, about 3 minutes
```

`run_sections_4_5.sh` and `run_section_7.sh` run the two halves separately and take the same option.
