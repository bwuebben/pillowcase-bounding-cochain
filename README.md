# The instanton homology of the (−2,3,q) pretzel knots and computed bounding cochains in the pillowcase

This repository contains two papers and the supporting computer code for results on both sides of
the **knot Atiyah–Floer program** for the pretzel family P(−2,3,q), q odd.

**Theorem (rigorous, unconditional).** rank I♮(P(−2,3,q)) = q + 2 for every odd q ≥ 3. The lower
bound is the Alexander-polynomial bound: the family's Alexander polynomials are Hironaka's
Lehmer-like polynomials (all coefficients in {0,±1}, exactly q+2 nonzero), derived here in closed
form by a validated Conway-skein recursion (`pillowcase/skein_alexander.py` — the recursion's unique
consistent normalization reproduces Lehmer's polynomial exactly at q=7). The upper bound is Manion's
closed-form reduced Khovanov homology of 3-strand pretzels (NYJM 24 (2018)): dim Khr = q + 2. The
Kronheimer–Mrowka squeeze closes the gap.

**Deficiency law (computed, perturbation-stable).** The naive (b = 0) pillowcase Lagrangian–Floer
homology of the family's natural tangle decomposition differs from I♮ by *exactly one differential*:

> naive rank − rank I♮ = 2·sgn(det K − 3),  det K = |q − 6|,

verified at q = 3, 5, 7, 11, 13 (naive ranks 5, 5, 7, 15, 17 vs I♮ = 5, 7, 9, 13, 15) at two
independent perturbations each. The direction flips at det = 3 — equivalently, where the double
branched cover stops being a homology sphere.

**The bounding cochains (computed; the first nonzero ones in pillowcase Floer theory).** The
correction conjectured by Cazassus–Herald–Kirk–Kotelskiy is realized in **both directions**:

> **P(−2,3,5)**: unique minimal b = s_A + s_B (two seam self-crossings); a μ³ **quadrilateral**
> cancels a bigon, rank 5 → 7.
> **P(−2,3,7)**: unique minimal b = a single self-crossing; a μ² **triangle** cancels a bigon,
> rank 7 → 9.
> **P(−2,3,11)**: fifty-five single-crossing cochains, each **creating** a differential,
> rank 15 → 13 — cancellation is rigid, creation is abundant.

## The two papers

- **`paper2/main.tex`** — *The instanton homology of the (−2,3,q) pretzel knots and computed
  bounding cochains in the pillowcase*. The main paper: the theorem, the deficiency law, and the
  three cochains. Figure 1 is generated from the computed curves (`pillowcase/make_figure.py` →
  `paper2/fig_q5.tex`).
- **`paper1/main.tex`** — *Trace-free SU(2) characters and ℤ/4 instanton gradings for
  two-bridge and (3,n)-torus knots* (companion). Assembles and independently verifies the
  representation-theoretic data (character varieties, ℤ/4 spectral-flow gradings, the first
  nonzero differential 8₁₉) underlying the pillowcase construction. Paper 2 takes up the
  bounding-cochain problem stated at the end of Paper 1.

Build either with `pdflatex main.tex` (run twice for cross-references); each compiles
independently with a standard TeX distribution. Compiled PDFs are included.

## The code

Everything is **pure Python 3** (standard library only — no NumPy, no dependencies). The
curves are built from first principles (quaternion representation theory), then validated
against Smith's published figures and numbers before the new computations are run. Each module
is also a self-checking test: run it directly and it prints a `PASS`/`FAIL` battery.

```bash
cd pillowcase

# --- validation gates: the reconstruction reproduces Smith's numbers (P(-2,3,5)) ---
python3 tangles.py     # conventions, Conway sum, seam fiber circles      (50/50)
python3 resolve.py     # the perturbed tangle sum / seam-circle resolution (11/11)
python3 earring.py     # the earring figure-eight; the 9-generator gate    (5/5)
python3 bigons.py      # winding-number bigon counter; 9 gens / 2 bigons / rank 5 (5/5)
python3 polygons.py    # generalized immersed-polygon counter              (2/2)

# --- the theorem's Alexander input ---
python3 skein_alexander.py   # sum|Delta(P(-2,3,q))| = q+2, Lehmer-validated   (13/13)

# --- the results ---
python3 b2_result.py         # P(-2,3,5): b = s_A + s_B, rank HF 5 -> 7 = I♮  (7/7, ~12s)
python3 pretzel_solve.py 3   # P(-2,3,7): b = one crossing (cancel), HF 7 -> 9 = I♮  (~3 min)
python3 pretzel_solve.py 5   # P(-2,3,11): 55 creation cochains, HF 15 -> 13 = I♮ (use do_quad=False for speed)
```

### Module guide

| file | role |
|---|---|
| `grounded.py`   | quaternion primitives (traceless SU(2) words) — the representation-theory base |
| `tangles.py`    | pillowcase coordinates, Conway sum as a fiber product, seam fiber circles |
| `resolve.py`    | Smith's cut-and-paste resolution of the seam circles → the **blue** Lagrangian |
| `earring.py`    | the Herald–Kirk earring figure-eight → the **red** Lagrangian; the 9-generator gate |
| `bigons.py`     | combinatorial (winding-number) bigon counter → the naive differential |
| `polygons.py`   | generalization to immersed (k+2)-gons (triangles, quadrilaterals, …) |
| `deform.py`     | the deformed differential ∂_b: triangle contributions; `build_pretzel(k)` for the family |
| `deform_full.py`| full ∂_b: triangles + quadrilaterals (with T²-consistent locality pruning) |
| `deform_pent.py`| the μ⁴ pentagon layer (found unnecessary — returns none) |
| `maurer_cartan.py`| the Maurer–Cartan data: μ⁰ (monogons), μ¹ (self-bigons), μ² (self-triangles) |
| `solve_b2.py`   | searches all supports for the rank-1 deformed differential; the unique minimal b (P(−2,3,5)) |
| `b2_result.py`  | **P(−2,3,5) result battery**: reproduces the gate, computes b, verifies rank 7 + Maurer–Cartan |
| `pretzel_solve.py`| **the family solver** `python3 pretzel_solve.py k`: builds P(−2,3,2k+1), searches for the MC-valid bounding cochain matching I♮ = q+2 (k=2: P(−2,3,5); k=3: P(−2,3,7); k=5: P(−2,3,11) creation direction) |
| `skein_alexander.py`| **Theorem input**: closed-form Alexander polynomials of the family via the validated Conway-skein Chebyshev recursion (Lehmer match at q=7); Σ\|Δ\| = q+2 |
| `make_figure.py`  | generates paper2's Figure 1 (`paper2/fig_q5.tex`) from the computed curves |
| `pert_check.py` | perturbation-stability: reruns the whole pipeline at a second perturbation |
| `diag_geometry.py`, `diag_cancel.py` | diagnostics used while developing the computation |

## Method and honesty

The computations are carried out **within the immersed-curve combinatorial model** of
Herald–Kirk and Smith, in which Lagrangian–Floer polygon counts are winding-number
computations for piecewise-linear immersed curves in the pillowcase. The values of `b` are the
values of the nonzero bounding cochains **under the standing conjecture** (of
Cazassus–Herald–Kirk–Kotelskiy) that this model computes the Floer-theoretic invariant — the
same hypothesis under which every computation in this subject is stated. The reconstruction is
validated against Smith's published output before the new computations run, and the results are
stable under change of perturbation.

The register is separated throughout: the **theorem** (rank I♮ = q+2) is proved unconditionally —
no model enters; the naive ranks, the deficiency law at its five members, and the cochains are
**computed within the model**; the law beyond q = 13 and the faithfulness of the model itself remain
**conjectural**. We claim no theorem about the analytic Atiyah–Floer correspondence. (The q = 17
member is untested: the current seam-resolution code cannot yet handle its overlapping seam arcs;
the law predicts naive rank 21.) See §1.5 and §6 of `paper2` for the precise statement.

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

## License

Code is released under the MIT License (`LICENSE`). The papers (`paper1/`, `paper2/`) are
© Bernd J. Wuebben; you may read and redistribute them for scholarly purposes with attribution.
