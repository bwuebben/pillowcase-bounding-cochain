# Shared programs: paper 2b, Section 3, and the checks of paper 1

Most programs in this directory carry out the computations of Section 3 of *Maurer–Cartan elements in the pillowcase and the Bar-Natan invariant for the (−2,3,q) pretzel knots* (`paper2b/`); the programs for its other
sections are in `paper2b/code/`. The last section below lists the numerical checks of paper 1.

For paper 2b, statement numbers refer to the compiled paper. They work in Smith's
coordinates with his perturbation for small t < 0 (Remark 3.2), unless the opposite sign is
named.

There are two kinds of computation, as in the paper.

* **Exact computations in the algebra 𝓑 of Kotelskiy, Watson and Zibrowius (Section 3.3, Propositions 3.7–3.12).** Each model
  curve is encoded as a reduced twisted complex over the algebra of Kotelskiy, Watson and
  Zibrowius. Morphism homology is computed in full by the mapping-cone reduction of Lemma 3.5,
  which raises an error unless every face complex has torsion homology. All arithmetic is over F₂.
* **Finite polygon tables (Appendix A.3, Computations A.1 and A.2).** These are the outputs of
  finite procedures and are recorded in the paper as data, not theorems. A matrix whose
  rank-derived statistic is displayed is first checked to square to zero; a matrix with nonzero
  square is not treated as a chain complex.

The code is pure Python 3 and uses only the standard library. The optional `--external` mode of
`q7_closure_probe.py` uses Spherogram 2.4.1 and Khoca 1.5.

## The exact q = 7 programs

| program | statements |
|---|---|
| `q7_kwz.py --encode` | Proposition 3.7: the encodings E (6 generators) and N (31 generators), the dimensions 7 and 9, the four switches b_s with D_N(b_s) = b_s² = 0; Corollary 3.9: the strict S18–S25 isomorphism |
| `q7_closure_probe.py --slope=-3/4 --selection-certificate` | Proposition 3.8: the second closure, its Alexander polynomial, the dimensions 31, 23, 25; Proposition 3.12: all sixteen sums in the span of the four switches are Maurer–Cartan, with their closure dimensions |
| `q7_quilt_census.py --all-node-census` | Proposition 3.10: the 82 self-intersection orbits (52 connector–main, 30 connector–connector); 73 generator-preserving census smoothings giving 45 distinct four-arrow Maurer–Cartan elements; the pair (9,31) only at S18 and S25 |
| `q7_quilt_census.py --two-switch-census --strict-pair-census` | Proposition 3.12: the 46 connector–main occurrences (38 distinct switches), the 703 sums of two switches (all Maurer–Cartan, 41 with pair (9,31)), none strictly isomorphic to the S18 presentation |
| `q7_excluded_orbits.py` | Proposition 3.11: the nine orbits where the census smoothing changes the module, with pairs (3,23), (5,23) three times and (7,25) five times; the opposite smoothings, with (9,31) only at S18 and S25, 21 generators and dimension 69 against the slope-3/4 earring, against 103 for (N, b_S18) |
| `c3_perturbed.py` | Proposition 3.4: the three transverse points of the partial fiber product over ((0,π/2),(0,π/2)), and the checks of the perturbed equations of Smith's C₃ configuration |
| `c3_q7_compare.py` | Remark 3.15: the numerical trace of Smith's perturbed curve at q = 7, its 50 self-intersection orbits, and the comparison of its encoding with the model |
| `surgery_check.py` | the comparison of the finite tables with single smoothings in Section 3 |

## The finite tables

| finite statement | program | output |
|---|---|---|
| q = 5 reconstruction: nine generators and two bigon entries | `tangles.py`, `resolve.py`, `earring.py`, `bigons.py` | the built-in checks pass |
| q = 5 bigon matrix | `b2_result.py` | D² = 0, rank 2, statistic 5 |
| q = 5 candidate support | `b2_result.py`, `solve_b2.py` | D² = 0, rank 1, statistic 7 |
| q = 7 singleton support | `pretzel_solve.py 3 --triangles-only --max-support 1` | D² = 0, statistic 9 |
| q = 11 singleton screen | `pretzel_solve.py 5 --triangles-only --max-support 1` | 62 candidates, 56 with D² = 0; the six failures are printed |
| q = 11, parameters (0.03, 0.10, 0.25) | `pretzel_solve.py 5 --triangles-only --max-support 1 --blue-epsilon .03 --red-epsilon .10 --red-pinch .25` | 66 candidates, 50 with D² = 0 |
| finite pentagon table | `deform_pent.py` | no contributions |
| two q = 5 perturbations | `pert_check.py` | as printed |
| closed form of the Alexander polynomials | `skein_alexander.py` | an independent check; the proof of Theorem 1.1 of paper 2 does not use it |

The two torus preimages of a self-intersection of a piecewise-linear curve are paired by the
involution with the tolerance `1e-7` (`maurer_cartan.orbit_group`). At q = 11 the model has 240
self-intersections, some of them less than `1e-3` apart.

## Module guide

| file | role |
|---|---|
| `grounded.py` | quaternion arithmetic for traceless SU(2) words |
| `tangles.py` | pillowcase coordinates, Conway sum, edge circles |
| `resolve.py` | the piecewise-linear resolution of the edge circles |
| `earring.py` | the piecewise-linear earring curve |
| `bigons.py`, `polygons.py` | finite bigon and polygon predicates |
| `deform.py`, `deform_full.py`, `deform_pent.py` | finite triangle, quadrilateral and pentagon tables; rank and square checks over F₂ |
| `maurer_cartan.py`, `solve_b2.py`, `pretzel_solve.py` | finite obstruction tables and support searches |
| `make_figure.py` | a TikZ drawing of the q = 5 curves, their generators and the two candidate-support crossings |

## Use

```sh
sh run_all.sh          # everything (the q = 11 screens take several minutes)
sh run_all.sh fast     # the reconstruction checks and the exact q = 7 programs
```

## Paper 1

The numerical checks of *Traceless SU(2) characters and ℤ/4 instanton gradings for two-bridge and
(3,n)-torus knots* (`paper1/`). No proof depends on them.

| program | checks |
|---|---|
| `riley_check.py` | the traceless Riley polynomials of the two-bridge knots with p ≤ 9, by exact computation of the Riley word |
| `torus_characters.py` | the count of traceless characters of T(3,n) for n ≤ 25, by enumerating arcs |
| `fs_gradings.py` | the ℤ/4 gradings of the irreducible characters of T(3,n) and the split between gradings 1 and 3 for odd n ≤ 43 |

