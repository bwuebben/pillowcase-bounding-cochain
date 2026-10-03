# Code for "Pillowcase Floer homology of the torus knots T(3,n)"

This directory contains the computations reported in the paper. Statement numbers refer to the
compiled paper.

There are two kinds of computation, as in the paper.

* **Rigorous interval arithmetic (Appendix A).** The three inequalities of Lemmas A.1–A.3 are proved
  by evaluating explicit functions on boxes in interval arithmetic (`mpmath.iv`, 30 decimal digits,
  π enclosed in an interval, outward rounding). Taylor coefficients are computed in exact rational
  arithmetic. These scripts are part of the proofs of Theorems 4.4 and 4.10 and Proposition 4.11.
* **Floating-point numerics (Sections 5–7).** The perturbed character varieties are traced by
  numerical continuation and the complexes are computed from the traced curves. These computations
  are numerical, not proofs. They confirm Theorem 6.10 in examples (Remark 6.12), and they are the
  basis of Proposition 7.1, which the paper labels as numerical.

Every script prints each number the paper states, next to the computed value, with `PASS` or
`FAIL`, and exits with status 0 if and only if all its checks pass.

## Scripts

| script | verifies | runtime |
|---|---|---|
| `interval_noncentral_circles.py` | Lemma A.1: F0 − 0.851\|F1\| ≥ 1.648 on Ω±, 96 × 96 = 9216 boxes, minimum 1.6482…; the independent Lipschitz check (240 × 240 centres, ≥ 1.637); the constant 1/(2 sin(π/5)) = 0.85065… of Theorem 4.4, attained at n = 7; the floating-point remark that F0 − c\|F1\| > 0 for every c < 2 | 25 s |
| `interval_central_circle.py` | Lemma A.2: the Taylor coefficients c3, c5, c7 (exact), the product-to-sum forms, the bounds M9 and \|R\|/ξ³, the lower bounds 22.06 and 19.18 for ξ ≤ 1/4, the 4096 boxes for ξ ≥ 1/4 with minima 0.798 and 0.650, and the independent check with splitting point ξ = 0.2 (22.80, 19.48) | 3 s |
| `interval_gamma_confinement.py` | Lemma A.3: N_j(ξ, 0) = 0, the factorization through Br_j (exact), the leading coefficients, the lower bounds 2.48 and 1.28 for ξ ≤ 0.3, and the adaptive boxes for ξ ≥ 0.3 (4100 and 4124 box evaluations, 1 and 7 boxes bisected, final covers of 4099 and 4117 boxes) | 1 s |
| `double_point_constants.py` | the inequality \|λ_c ε_A\| > 2\|μ_c ε_B\| of Remark 5.8 at the sampled perturbations (n ≡ 5 mod 6, n ≤ 47); the values after Lemma 4.13 (V_γ = 2/√3, V_φ = 2√3 for n = 4; formula against direct evaluation for n = 4, 5, 10, 11, 16, 17, 41) and after Lemma 5.4 (μ_c, λ_c for n = 4, 5, 10, 11, 16, 17, 22, 23) | < 1 s |
| `perturbed_complexes.py` | Remark 6.12 and Remark 5.8: the twelve values n = 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 22, 23 at the two perturbations (nε_A, ε_B, ε) = (0.3, 0.05, 0.01), (−0.3, 0.03, 0.02): arc endpoint and lines crossed (Table 1, Lemma 6.1), generators (Proposition 6.3), bigons (Proposition 6.5), degrees from formula (2.2) (Proposition 6.9), graded homology (Table 2), the comparison with Hedden–Herald–Kirk [HHK18, §§11.2, 11.6, 11.7] for n = 4, 5, 7, and the tangency of L1 with the line field of slope one; the resolution observed at these perturbations (split for n ≡ 4, joined for n ≡ 5 mod 6); and, for Remark 5.8, the same checks at three perturbations on the ε_B-axis (n = 5 at (ε_A, ε_B) = (0, 0.04), n = 11, 17 at (0, 0.03)), where the resolution is split for n ≡ 5 | 75 s |
| `survey_n_le_50.py` | Remark 6.12, last part: generator counts, bigons and homology ranks for all 4 ≤ n ≤ 50 coprime to 3 at four perturbations each (128 cases); the resolution observed at these perturbations (Remark 6.12) | 6 min |
| `t34_second_decomposition.py` | Section 7 and Proposition 7.1: T(3,4) with (r,s) = (7,−5) at 36 perturbations (twelve directions, sizes 0.03, 0.01, 0.003), items (1)–(5); the polygon representative of the arc; the twelve self-intersections; the bounding cochain s and its triangle; the controls (r,s) = (3,−2) and (−1,1); the strand-meridian quotients. The approximate coordinates quoted in Section 7 (the crossings at γ ≈ 4.33π and 2.67π, the point s ≈ (0.18π, 0.10π)) are checked over all 36 runs; an additional run at (ε_A, ε_B) = (0.075, 0.05) checks the triangle and the self-intersection count | 4.5 min |

Shared code:

* `t3n_pillowcase.py`: the equations of the tangle of Hedden, Herald and Kirk [HHK18, (31)–(34)], the restriction map to the
  pillowcase, continuation, lifts to the plane, the earring curve and its strands, generators,
  candidate bigons (Lemma 6.4), and the grading formula (2.2).
* `trig_series.py`: exact arithmetic for trigonometric polynomials whose frequencies are linear in h,
  used for the Taylor coefficients in Lemmas A.2 and A.3.

## Requirements and use

Python 3.11 or later (tested with 3.13), with the packages in `requirements.txt` (numpy and
mpmath 1.3.0). For example:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PY=.venv/bin/python ./run_all.sh
```

Each script can also be run on its own, e.g. `.venv/bin/python interval_central_circle.py`.
`./run_all.sh --quick` runs a shorter version of the two longest computations (the Section 7 sweep
at size 0.03 only, and the survey up to n = 20). `perturbed_complexes.py` accepts a list of values
of n.

Runtimes above are wall-clock times on a laptop (Apple silicon, one core). The total for
`./run_all.sh` is about 12 minutes.
