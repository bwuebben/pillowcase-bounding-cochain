# Computed bounding cochains in the pillowcase: P(−2,3,5) and P(−2,3,7)

This repository contains two papers and the supporting computer code for explicit
computations on the **pillowcase (symplectic) side of the knot Atiyah–Floer conjecture**.

The headline results are the **first explicitly computed nonzero pillowcase bounding
cochains**. The naive Lagrangian–Floer homology of a knot's tangle decomposition is not a
knot invariant; the conjectural repair of Cazassus–Herald–Kirk–Kotelskiy is a *bounding
cochain* `b` supported on the self-intersections of the immersed curves, satisfying a mod-2
Maurer–Cartan equation. Kai Smith ([arXiv:2412.06066](https://arxiv.org/abs/2412.06066))
proved that such a bounding cochain must be nonzero for P(−2,3,5), but computed no value, and
no nonzero pillowcase bounding cochain had been computed. This code computes two of them —
the two smallest members of the pretzel family P(−2,3,2k+1):

> **P(−2,3,5) = T(3,5).** `b = s_A + s_B`, the sum of two self-intersection points of the
> perturbed Conway-sum Lagrangian 𝓡ₜ(Q₁ᐟ₃ + Q₁ᐟ₅), one on each pillowcase seam (at
> s_A ≈ (γ,θ) = (0.03, 1.27) and s_B ≈ (3.06, 4.98)). It raises the rank from the naive **5**
> to **7 = rank I♮**, via a single immersed **quadrilateral** (μ³) that cancels one bigon.
>
> **P(−2,3,7).** `b` is a **single** self-intersection at ≈ (0.05, 5.41). It raises the rank
> from the naive **7** to **9 = rank I♮**, via a single immersed **triangle** (μ²) that
> cancels one bigon.

In each case `b` is the unique minimal bounding cochain matching I♮, and cancels exactly one
bigon (+2 to the homology rank). The two examples realize the **same +2 correction through
different polygon orders** (quadrilateral vs. triangle) — the paper's central structural
observation. Both values are stable under change of perturbation.

## The two papers

- **`paper2/main.tex`** — *Computed bounding cochains in the pillowcase: the pretzel knots
  P(−2,3,5) and P(−2,3,7)*. The main results. Reconstructs Smith's Lagrangians from the
  quaternionic representation theory, reproduces his 9 generators / 2 bigons / rank-5 gate as
  validation, then computes the two bounding cochains, and discusses the family honestly.
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

# --- the results ---
python3 b2_result.py     # P(-2,3,5): b = s_A + s_B, rank HF 5 -> 7 = I♮   (7/7, ~12s)
python3 pretzel_solve.py 3   # P(-2,3,7): b = one crossing, rank HF 7 -> 9 = I♮  (~3 min)
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
| `pretzel_solve.py`| **the family solver** `python3 pretzel_solve.py k`: builds P(−2,3,2k+1), searches for the MC-valid bounding cochain matching I♮ (k=2 reproduces P(−2,3,5); k=3 gives P(−2,3,7)) |
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

We claim **no theorem** about the analytic Atiyah–Floer correspondence, and **no uniform
statement** about the family P(−2,3,2k+1) beyond its two smallest members: the computed data
(paper §5.3) shows the naive-rank pattern does not obviously persist (P(−2,3,11) has
perturbation-stable naive rank 15), and rank I♮ for q ≥ 11 is not determined here. See §1.4 and
§5 of `paper2` for the precise separation of what is proved, what is computed within the model,
and what remains conjectural.

## References

- K. Smith, *Perturbed traceless SU(2) character varieties of tangle sums*, arXiv:2412.06066.
- G. Cazassus, C. Herald, P. Kirk, A. Kotelskiy, *The correspondence induced on the pillowcase
  by the earring tangle*, J. Topol. 15 (2022), arXiv:2010.04320.
- M. Hedden, C. Herald, P. Kirk, *The pillowcase and traceless representations of knot groups
  I, II*, arXiv:1301.0164, arXiv:1501.00028.
- C. Herald, P. Kirk, *An endomorphism on immersed curves in the pillowcase*, arXiv:2407.11247.
- M. Akaho, D. Joyce, *Immersed Lagrangian Floer theory*, J. Differential Geom. 86 (2010),
  arXiv:0803.0717.

## License

Code is released under the MIT License (`LICENSE`). The papers (`paper1/`, `paper2/`) are
© Bernd J. Wuebben; you may read and redistribute them for scholarly purposes with attribution.
