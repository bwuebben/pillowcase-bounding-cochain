# A computed bounding cochain: the pretzel knot P(−2,3,5)

This repository contains two papers and the supporting computer code for an explicit
computation in the **pillowcase (symplectic) side of the knot Atiyah–Floer conjecture**.

The headline result is the first explicitly computed **nonzero pillowcase bounding
cochain**. For the pretzel knot P(−2,3,5) = T(3,5), Kai Smith
([arXiv:2412.06066](https://arxiv.org/abs/2412.06066)) proved that the bounding cochain of
one tangle in a natural decomposition must be nonzero — by showing the naive Lagrangian–Floer
homology has rank 5 while the reduced singular instanton homology I♮ has rank 7 — but he did
not compute it. This code computes it:

> **b₂ = s_A + s_B**, the sum of two self-intersection points of the perturbed Conway-sum
> Lagrangian 𝓡ₜ(Q₁ᐟ₃ + Q₁ᐟ₅), one on each pillowcase seam (at coordinates
> s_A ≈ (γ,θ) = (0.03, 1.27) and s_B ≈ (3.06, 4.98)).
>
> It is the **unique minimal** bounding cochain for which the deformed Floer homology has
> rank **7 = rank I♮(P(−2,3,5))**, raising the naive rank of 5 to Kronheimer–Mrowka's 7.

The mechanism is a single immersed **quadrilateral** (a higher A∞ product μ³) whose two
interior vertices are s_A, s_B and which cancels one of the two naive bigons modulo 2. The
Maurer–Cartan equation is satisfied automatically because the blue Lagrangian's self-Floer
products μ⁰, μ¹, μ² all vanish.

## The two papers

- **`paper2/main.tex`** — *A computed bounding cochain: the pretzel knot P(−2,3,5)*. The main
  result. Reconstructs Smith's two Lagrangians from the quaternionic representation theory,
  reproduces his 9 generators / 2 bigons / rank-5 gate as validation, then computes b₂.
- **`paper1/main.tex`** — *Trace-free SU(2) characters and ℤ/4 instanton gradings for
  two-bridge and (3,n)-torus knots* (companion paper). Assembles and independently verifies
  the representation-theoretic data (character varieties, ℤ/4 spectral-flow gradings, the
  first nonzero differential 8₁₉) underlying the pillowcase construction. Paper 2 resolves the
  open problem stated at the end of Paper 1.

Build either with `pdflatex main.tex` (run twice for cross-references); each compiles
independently and needs only a standard TeX distribution.

## The code

Everything is **pure Python 3** (standard library only — no NumPy, no dependencies). The
curves are built from first principles (quaternion representation theory), then validated
against Smith's published figures and numbers before the new computation is run. Each module
is also a self-checking test: run it directly and it prints a `PASS`/`FAIL` battery.

```bash
cd pillowcase

# --- validation gates: the reconstruction reproduces Smith's numbers ---
python3 tangles.py     # conventions, Conway sum, seam fiber circles      (50/50)
python3 resolve.py     # the perturbed tangle sum / seam-circle resolution (11/11)
python3 earring.py     # the earring figure-eight; the 9-generator gate    (5/5)
python3 bigons.py      # winding-number bigon counter; 9 gens / 2 bigons / rank 5 (5/5)
python3 polygons.py    # generalized immersed-polygon counter              (2/2)

# --- the result ---
python3 b2_result.py   # THE RESULT: b₂ = s_A + s_B, rank HF = 7 = I♮      (7/7, ~12s)
```

### Module guide

| file | role |
|---|---|
| `grounded.py`   | quaternion primitives (traceless SU(2) words) — the representation-theory base |
| `tangles.py`    | pillowcase coordinates, Conway sum as a fiber product, seam fiber circles |
| `resolve.py`    | Smith's cut-and-paste resolution of the seam circles → the **blue** Lagrangian |
| `earring.py`    | the Herald–Kirk earring figure-eight → the **red** Lagrangian; the 9-generator gate |
| `bigons.py`     | combinatorial (winding-number) bigon counter → the naive rank-5 differential |
| `polygons.py`   | generalization to immersed (k+2)-gons (triangles, quadrilaterals, …) |
| `deform.py`     | the deformed differential ∂_b: triangle contributions through each self-crossing |
| `deform_full.py`| full ∂_b: triangles + quadrilaterals (with T²-consistent locality pruning) |
| `deform_pent.py`| the μ⁴ pentagon layer (found unnecessary — returns none) |
| `maurer_cartan.py`| the Maurer–Cartan obstruction: μ⁰ (monogons), μ¹ (self-bigons), μ² (self-triangles) — all vanish |
| `solve_b2.py`   | searches all supports for rank(∂_b) = 1; finds the unique minimal b₂ |
| `b2_result.py`  | **end-to-end result battery**: reproduces the naive gate, computes b₂, verifies rank 7 and Maurer–Cartan |
| `pert_check.py` | perturbation-stability: reruns the whole pipeline at a second perturbation |
| `diag_geometry.py`, `diag_cancel.py` | diagnostics used while developing the computation |

## Method and honesty

The computation is carried out **within the immersed-curve combinatorial model** of
Herald–Kirk and Smith, in which Lagrangian–Floer polygon counts are winding-number
computations for piecewise-linear immersed curves in the pillowcase. The value of b₂ is the
value of Smith's proven-nonzero bounding cochain **under the standing conjecture** (of
Cazassus–Herald–Kirk–Kotelskiy) that this model computes the Floer-theoretic invariant — the
same hypothesis under which every computation in this subject is stated. The reconstruction of
the two Lagrangians is validated against Smith's published output before the new computation
is run, and the result is stable under change of perturbation. See §1.4 and §5 of `paper2`
for the precise separation of what is proved, what is computed within the model, and what
remains conjectural.

## References

- K. Smith, *Perturbed traceless SU(2) character varieties of tangle sums*, arXiv:2412.06066.
- G. Cazassus, C. Herald, P. Kirk, A. Kotelskiy, *The correspondence induced on the pillowcase
  by the earring tangle*, J. Topol. 15 (2022), arXiv:2010.04320.
- M. Hedden, C. Herald, P. Kirk, *The pillowcase and traceless representations of knot groups
  I, II*, arXiv:1301.0164, arXiv:1501.00028.
- C. Herald, P. Kirk, *An endomorphism on immersed curves in the pillowcase*, arXiv:2407.11247.

## License

Code is released under the MIT License (`LICENSE`). The papers (`paper1/`, `paper2/`) are
© Bernd J. Wuebben; you may read and redistribute them for scholarly purposes with attribution.
