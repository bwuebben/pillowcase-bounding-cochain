# Atiyah–Floer pillowcase papers and exact computations

This repository contains three papers and the supporting computer code for results on both sides of
the **knot Atiyah–Floer program** for the pretzel family P(−2,3,q), q odd.

**Theorem (rigorous, unconditional, integral).** For every odd q ≥ 3, the reduced singular
instanton knot homology I♮(P(−2,3,q)) is free abelian of rank q + 2. The proof combines Hironaka's
exact formula for the family's Alexander polynomials (Lehmer-like: all coefficients in {0,±1},
exactly q+2 nonzero), Manion's closed-form reduced integral Khovanov homology of 3-strand pretzels
(NYJM 24 (2018)), and the integral Kronheimer–Mrowka spectral sequence. The Alexander input is
re-derived here independently by a validated Conway-skein recursion
(`pillowcase/skein_alexander.py` — its unique consistent normalization reproduces Lehmer's
polynomial exactly at q = 7).

**The q = 7 Maurer–Cartan calculation.** In the wrapped Fukaya subcategory that CHKK identify
with twisted complexes over the two-arc algebra of Kotelskiy–Watson–Zibrowius, the higher products
vanish, so the Maurer–Cartan equation for a deformation is the finite identity (δ+b)² = 0. The
morphism complex over that algebra is infinite dimensional; a filtration lemma (Paper 2, Lemma 4.1)
reduces its homology to three integers for any pair of finite twisted complexes, without truncating
word length. The resolved `Q_{1/3}+Q_{1/7}` curve encodes to a 31-generator twisted complex; the
pairing with the earring has dimension 7, and smoothing the curve at any one of four
self-intersection orbits (`S18`, `S25`, `S69`, `S74`) gives a four-arrow Maurer–Cartan element
raising it to 9. The four objects represent exactly three homotopy classes, with `S18` and `S25`
strictly isomorphic. A second closure of the same tangle, `num(Q_{-3/4}+Q_{1/3}+Q_{1/7})`, has
instanton rank 31 and pairing dimensions 31, 23 and 25 against the three classes. An exhaustive
enumeration of all 82 self-intersection orbits — 52 connector–main and 30 connector–connector —
finds 73 generator-preserving smoothings determining 45 distinct, linearly independent four-arrow
elements; the two closure ranks select only `S18/S25` among them. Within the connector–main
subfamily (46 occurrences, 38 distinct elements) every element of the `2^38`-element span is
Maurer–Cartan, and 41 of the `C(38,2) = 703` two-element sums share the same two ranks, so the
selection genuinely depends on a local-support hypothesis. Separately, the Lagrangian
correspondence induced by the `Q_{1/3}` line is immersed with a triple point, hence not embedded,
so Gao's wrapped representability theorem does not apply to it.

**Paper III exact-action theorem.** For the pillowcase composition associated to
`Q_{1/3}+Q_{1/7}`, the small-perturbation correspondence is exact on generic symmetric compact
truncations, and a half-period symmetry fixes the relative primitive between the two output arcs.
The singular path word has 58 seam-overlap candidates, exactly 50 of which persist for small
positive perturbation: 28 diagonal–circle and 22 circle–circle roots. Their limiting actions are
computed exactly as rational multiples of π². Forty-two are negative, eight are positive, none
vanishes, and the minimum absolute limiting action is π²/42. The complete census and all signs
persist for sufficiently small positive perturbation. This theorem makes no holomorphic-quilt
count, compactification, bounding-cochain, or instanton–pillowcase comparison claim.

**Finite pillowcase computations.** The pillowcase code reports finite polygon computations for
explicitly specified piecewise-linear immersed curves: finite bigon matrices, triangle and
quadrilateral tables inside configured edge windows, and screens for candidate correction supports.
The finite bigon statistic h = n − 2·rank(D) takes the values 5, 5, 7, 15, 17 at
q = 3, 5, 7, 11, 13 (the q = 3 value is Hedden–Herald–Kirk's; this pipeline needs gcd(3,q) = 1),
against the instanton rank q + 2 = 5, 7, 9, 13, 15 — a difference of at most two, whose sign flips
at det K = |q − 6| = 3. The low-order tables single out candidate supports whose finite matrices
square to zero over GF(2) and have the target statistic at the reported perturbation: a
two-crossing support at q = 5 (h = 7), a one-crossing support at q = 7 (h = 9), and, at q = 11,
fifty-five singleton supports passing the finite screen, of which fifty-two square to zero
(56 and 42 at an independent perturbation).

**What is not claimed.** The q = 5, 11, 13 candidate tables remain statistics of explicitly
constructed finite matrices; their supports are not proved to satisfy the full immersed-Fukaya
Maurer–Cartan equation, and the finite edge windows are not proved exhaustive. The q = 7
Maurer–Cartan statements do avoid those algebraic truncations, but they are statements about
explicitly presented objects in `Tw(B)`. They do not prove that the CHKK instanton–pillowcase
tangle assignment exists or selects one of the stated classes. Paper 2 isolates that remaining
local-support statement as a numbered Hypothesis.

**Revision of August 2026 (arXiv v3).** Earlier versions of the papers (arXiv v1–v2) and of this README
described these outputs as computed bounding cochains with unique minimal supports and as deformed
Floer homology values. Those claims are withdrawn. The earlier q = 11 screen did not test
square-zero — three of its fifty-five matrices fail it — and the finite searches were never proved
exhaustive. The v3 Paper 2 adds the q = 7 Maurer–Cartan calculation with its finiteness lemma, Smith's
regular-homotopy identification of the underlying main immersion, the exhaustive enumeration of
single smoothings, the triple-point obstruction to Gao's representability theorem, and a numbered
statement of the remaining local-support hypothesis. It does not reinstate the withdrawn
finite-polygon claims. Both papers carry an Appendix A itemizing every change from v2.

## The three papers

- **`paper3/main.tex`** (**full draft, August 15, 2026**) — *Exact action separation for a
  pretzel-tangle composition in the pillowcase*. Proves exactness and half-period normalization for
  the `Q_{1/3}+Q_{1/7}` composition, the six-fold and 21-sheet description, the exact 50-root
  census, the 42-negative/eight-positive action split, and small-perturbation persistence. It
  explicitly stops before any quilt count or bounding-cochain conclusion.
- **`paper2/main.tex`** (**arXiv:2607.26096**) — *The instanton homology of the (−2,3,q) pretzel
  knots and Maurer–Cartan deformations in the two-arc algebra of the pillowcase*. The main paper:
  the integral instanton theorem, the finiteness lemma for morphism complexes over the two-arc
  algebra, the q = 7 Maurer–Cartan and second-closure calculations, Smith's identification of the
  underlying main immersion, the triple-point obstruction, the finite candidate-support tables for
  comparison, and the precise remaining local-support hypothesis. Figure 1 is generated from the computed curves
  (`pillowcase/make_figure.py`).
- **`paper1/main.tex`** (**arXiv:2607.26095**) — *Traceless SU(2) characters and ℤ/4 instanton
  gradings for two-bridge and (3,n)-torus knots* (companion). Determines the ℤ/4 gradings of the
  (3,n)-torus knots through the double branched cover rather than an index computation — the two
  traceless lifts of a flat connection carry the same grading, so Daemi–Scaduto's torus-knot theorem
  transports to a grading split — recovering the chain-rank distribution conjectured by
  Poudel–Saveliev and established for all torus knots by Daemi–Scaduto; and comparing that rank
  vector with the Alexander norm gives the total rank of the framed differential in every residue
  class of n mod 6 (zero for n ≡ 1, 2, one for n ≡ 4, 5). Also the two-bridge traceless characters
  via the Riley polynomial, the (3,n) dihedral dichotomy, and the 8₁₉ differential. Paper 2 takes up the
  correction problem stated at the end of Paper 1.

Build any paper with `pdflatex main.tex` (run twice for cross-references); each compiles
independently with a standard TeX distribution. Compiled PDFs are included under their
distinctive names (`traceless-gradings.pdf`, `pretzel-cochains.pdf`,
`exact-action-separation.pdf`).

## The methods note

- **`note/main.tex`** (compiled: `note/cochain-search-note.pdf`) — *Finding the pillowcase
  bounding cochains: the search behind the P(−2,3,q) computations, in complete detail* (12 pp).
  A correspondence note, written for Zuyi Zhang, expanding §§4–5 of paper 2 as of July 2026 to
  full operational detail: the lift bookkeeping, the pruning discipline, the traps that produced
  wrong answers along the way, and the acceptance criteria, together with what can be said about
  the analytic question.

  **Status (August 12, 2026).** The note predates the August 12 revision of paper 2 and presents
  the search as it was then understood. In particular its statements that the q = 11 screen
  produced fifty-five cochains and that the searches establish minimality by exhaustion are
  superseded — see the revision note above and the current paper 2. It is kept unchanged as part
  of the correspondence record.

## The code

Everything is **pure Python 3** (standard library only — no NumPy, no dependencies). The
curves are built from first principles (quaternion representation theory), then validated
against Smith's published figures and numbers before the new computations are run. Each module
is also a self-checking test: run it directly and it prints a `PASS`/`FAIL` battery. Every matrix
whose rank statistic is displayed is first subjected to a GF(2) square-zero audit; a matrix with
nonzero square is not treated as a chain complex.

```bash
cd pillowcase

# --- validation gates: the reconstruction reproduces Smith's numbers (P(-2,3,5)) ---
python3 tangles.py     # conventions, Conway sum, seam fiber circles      (50/50)
python3 resolve.py     # the perturbed tangle sum / seam-circle resolution (11/11)
python3 earring.py     # the earring figure-eight; the 9-generator gate    (5/5)
python3 bigons.py      # winding-number bigon counter; 9 gens / 2 bigons   (5/5)
python3 polygons.py    # generalized immersed-polygon counter              (2/2)

# --- the theorem's Alexander input (independent regression; the proof cites Hironaka) ---
python3 skein_alexander.py   # sum|Delta(P(-2,3,q))| = q+2, Lehmer-validated   (13/13)

# --- paper 1's verifications ---
python3 riley_check.py       # Thm 1.1: the traceless Riley polynomial, exactly over Z[i][u]  (91 checks)
python3 torus_characters.py  # Prop 4.1 / Thm 1.2: the (3,n) traceless character count        (65 checks)
python3 fs_gradings.py       # Z/4 spectral-flow gradings of IC♮(T(3,n)), Anvari-verified

# --- q=7 Maurer–Cartan and second-closure computations ---
python3 q7_kwz.py --encode
python3 q7_closure_probe.py --slope=-3/4 --selection-certificate
python3 q7_quilt_census.py --pairing-census
python3 q7_quilt_census.py --two-switch-census --strict-pair-census

# --- Paper III: exact singular-limit root and action certificate ---
python3 q7_exact_actions.py

# --- actual corrected-C3 numerical stability diagnostic (not an analytic proof) ---
python3 c3_perturbed.py
python3 c3_q7_compare.py --stability --switch-census

# --- the finite candidate-support computations ---
python3 b2_result.py         # q=5: support {s_A,s_B}; D^2=0, h: 5 -> 7   (~15 s)
python3 pretzel_solve.py 3 --triangles-only --max-support 1
                             # q=7: default singleton; D^2=0, h: 7 -> 9
python3 pretzel_solve.py 5 --triangles-only --max-support 1
                             # q=11: 55 pass the finite screen; 52 square to zero (3 failures printed)
python3 maurer_cartan.py     # finite monogon/self-bigon/self-triangle tables (diagnostics only)
python3 pert_check.py        # two finite q=5 perturbation runs
python3 surgery_check.py     # smoothing blue at the q=5,7 supports and recounting bigons
                             # reproduces the finite deformed matrix entrywise (~2 min)
```

### Module guide

| file | role |
|---|---|
| `grounded.py`   | quaternion primitives (traceless SU(2) words) — the representation-theory base |
| `tangles.py`    | pillowcase coordinates, Conway sum as a fiber product, seam fiber circles |
| `resolve.py`    | Smith's cut-and-paste resolution of the seam circles → the **blue** curve |
| `earring.py`    | the Herald–Kirk earring figure-eight → the **red** curve; the 9-generator gate |
| `bigons.py`     | finite (winding-number) bigon predicate and matrix → the undeformed tables |
| `polygons.py`   | finite immersed (k+2)-gon predicate (triangles, quadrilaterals, …) |
| `deform.py`     | triangle tables, curve construction, reusable GF(2) rank and square-zero checks; `build_pretzel(k)` for the family |
| `deform_full.py`| triangle + distinct-pair quadrilateral tables (finite windows) |
| `deform_pent.py`| distinct-support finite pentagon table |
| `maurer_cartan.py`| finite obstruction-table diagnostics: monogons, self-bigons, self-triangles |
| `solve_b2.py`   | finite support search over the cached tables (q = 5) |
| `b2_result.py`  | self-checking q = 5 candidate battery: the two-crossing support, D²=0, h = 7 |
| `pretzel_solve.py`| the family candidate solver `python3 pretzel_solve.py k` with a mandatory GF(2) square-zero audit (k=2: q=5; k=3: q=7; k=5: q=11) |
| `q7_kwz.py` | the 31-generator q = 7 twisted complex over the two-arc algebra, the four Maurer–Cartan elements from local smoothings, the finite reduction of the wrapped morphism complex, the strict `S18`–`S25` isomorphism, and the three component classes |
| `q7_closure_probe.py` | slope `-3/4` rational-earring pairing, auxiliary Alexander/Khovanov certificate, and the complete sixteen-element named-switch span |
| `q7_quilt_census.py` | exhaustive enumeration of the 82 self-intersection orbits (73 generator-preserving smoothings, 45 distinct elements) and its 52-orbit connector–main subfamily (46 occurrences, 38 elements), all 703 two-element sums, closure-rank distribution, and strict relabeling test |
| `q7_exact_actions.py` | Paper III's standard-library rational certificate: singular path word, seven algebraic normal-root boxes, 58 overlap candidates, 50 persistent roots, the 25-row area–residue table, and the exact 42/8 sign split |
| `c3_perturbed.py` | corrected equation-level trace of Smith's `C3` correspondence; numerical diagnostic with explicit source caveats |
| `c3_q7_compare.py` | two-perturbation comparison of the actual `C3` trace with the PL twisted complex and its 43-singleton enumeration |
| `skein_alexander.py`| closed-form Alexander polynomials of the family via the validated Conway-skein Chebyshev recursion (Lehmer match at q=7); Σ\|Δ\| = q+2 — an independent regression; the theorem cites Hironaka |
| `make_figure.py`  | generates paper2's Figure 1 from the computed curves |
| `riley_check.py`  | **paper 1, Thm 1.1**: the Riley word of 𝔟(p,q) at meridian eigenvalue s = i, computed exactly in Z[i][u]; the traceless Riley polynomial as a gcd over Q(i), matched against φ_p = ∏(u + 4sin²(πk/p)) built as an exact integer polynomial. Also proves the binary-dihedral identity A(AB)A⁻¹ = (AB)⁻¹ symbolically |
| `torus_characters.py` | **paper 1, Prop 4.1 / Thm 1.2**: direct enumeration of the (3,n) representation arcs and their traceless characters for all n ≤ 25 (independent check of the closed-form count, proved in the paper); the dihedral dichotomy; and the cross-check N(3,n) = 2a against `fs_gradings.py` for all odd n ≤ 43 |
| `fs_gradings.py`  | **paper 1**: Fintushel–Stern / equivariant-ρ spectral-flow gradings assembling the framed chain complex of IC♮(T(3,n)), verified against Anvari Ex. 6.1 |
| `pert_check.py` | perturbation stability of the finite tables: reruns the q = 5 pipeline at a second perturbation |
| `surgery_check.py` | machine check of the combinatorial surgery statement at q = 5, 7: smoothing blue ι-equivariantly at the support crossings and recounting bigons reproduces the finite deformed matrix entrywise; the other sector choice reproduces the undeformed matrix |
| `diag_geometry.py`, `diag_cancel.py` | diagnostics used while developing the computation |

## Method and scope of claims

The pillowcase computations are carried out in the piecewise-linear immersed-curve model of
Herald–Kirk and Smith, with polygon counts implemented as winding-number computations inside
configured edge windows and up to a configured polygon order. The reconstruction is validated
against Smith's published output before the new computations run.

The scope is separated throughout. The **instanton theorem** (I♮ free abelian of rank q + 2) is
proved unconditionally from cited results — no pillowcase model enters. The q = 7 Maurer–Cartan and
second-closure propositions are algebraic statements in the specified two-arc category, verified
without truncation of word length.
Paper III is likewise unconditional but narrower: it proves the finite exact geometry and action
filtration of one composed curve. Its final energy statement is only an obstruction for any exact
configuration already satisfying Stokes' identity; it does not assert that a quilted moduli space
exists or that any positive-action root is counted.
The remaining identification of that selected deformation with the instanton tangle object is
conditional on the CHKK correspondence and the q = 7 figure-eight local-support hypothesis. The
other bigon statistics and candidate supports are finite computations within the model,
square-zero checked over GF(2), not claimed Floer differentials. We claim no theorem about the
analytic Atiyah–Floer correspondence. (The q = 17 member is untested: the current seam-resolution
code cannot yet handle its overlapping seam arcs.)

## References

- K. Smith, *Perturbed traceless SU(2) character varieties of tangle sums*, arXiv:2412.06066.
- G. Cazassus, C. Herald, P. Kirk, A. Kotelskiy, *The correspondence induced on the pillowcase
  by the earring tangle*, J. Topol. 15 (2022), arXiv:2010.04320.
- M. Hedden, C. Herald, P. Kirk, *The pillowcase and traceless representations of knot groups
  I, II*, arXiv:1301.0164, arXiv:1501.00028.
- C. Herald, P. Kirk, *An endomorphism on immersed curves in the pillowcase*, arXiv:2407.11247.
- M. Akaho, D. Joyce, *Immersed Lagrangian Floer theory*, J. Differential Geom. 86 (2010),
  arXiv:0803.0717.
- A. Manion, *The Khovanov homology of 3-strand pretzels, revisited*, New York J. Math. 24 (2018),
  1076–1100, arXiv:1303.3303.
- E. Hironaka, *The Lehmer polynomial and pretzel links*, Canad. Math. Bull. 44 (2001), 440–451.
- Y. Lim, *Instanton homology and the Alexander polynomial*, Proc. Amer. Math. Soc. 138 (2010),
  3759–3768.
- A. Daemi, C. Scaduto, *Chern–Simons functional, singular instantons, and the four-dimensional
  clasp number*, J. Eur. Math. Soc. 26 (2024), 2127–2190, arXiv:2007.13160.
- P. Poudel, N. Saveliev, *Link homology and equivariant gauge theory*, Algebr. Geom. Topol. 17
  (2017), 2635–2687, arXiv:1502.03116.
- A. Daemi, C. Scaduto, *Equivariant aspects of singular instanton Floer homology*, Geom. Topol.
  28 (2024), 4057–4190, arXiv:1912.08982.
- A. Kotelskiy, L. Watson, C. Zibrowius, *Immersed curves in Khovanov homology*, arXiv:1910.14584.
- Y. Gao, *Functors of wrapped Fukaya categories from Lagrangian correspondences*, arXiv:1712.00225.
- D. Schütz, *On the Khovanov homology of 3-braids*, Quantum Topol., published online November 7,
  2025, doi:10.4171/QT/248, arXiv:2501.11547.

## License

Code is released under the MIT License (`LICENSE`). The papers and the note (`paper1/`,
`paper2/`, `paper3/`, `note/`) are
© Bernd J. Wuebben; you may read and redistribute them for scholarly purposes with attribution.
