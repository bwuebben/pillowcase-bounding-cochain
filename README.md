# Atiyah–Floer pillowcase papers and exact computations

This repository contains four papers on the knot version of the **Atiyah–Floer correspondence**,
which relates Kronheimer–Mrowka's reduced singular instanton homology I♮(K) to Lagrangian Floer
homology of immersed curves in the pillowcase, the traceless SU(2) character variety of the
four-punctured sphere. It also contains the computer code behind the computations. Two families
are treated: the **torus knots** (papers 1, 4 and 5) and the **pretzel knots P(−2,3,q)** (paper 2).
Each paper cites the others where their results meet.

| directory | paper | pages |
|---|---|---|
| [`paper1/`](paper1) | *Traceless SU(2) characters and ℤ/4 instanton gradings for two-bridge and (3,n)-torus knots* (an earlier version is arXiv:2607.26095) | 14 |
| [`paper4/`](paper4) | *The S-complex differential of torus knots* | 12 |
| [`paper5/`](paper5) | *Pillowcase Floer homology of the torus knots T(3,n)* | 40 |
| [`paper2/`](paper2) | *Instanton and pillowcase homology of the (−2,3,q) pretzel knots* (an earlier version, under a different title, is arXiv:2607.26096) | 60 |

For the torus knots, read papers 1, 4 and 5 in that order. Paper 1 assembles the instanton-side
facts, paper 4 proves new results about the S-complex differential for every torus knot, and paper 5
computes the pillowcase side for the whole T(3,n) family. Paper 2 treats the pretzel family; its
Sections 4 and 5 build on the curves of paper 5, and its first two members, P(−2,3,3) = T(3,4) and
P(−2,3,5) = T(3,5), are torus knots covered by papers 1, 4 and 5.

## The papers

### Paper 1: traceless characters and instanton gradings (`paper1/`)

An account of the generators, the ℤ/4 gradings and the differential of Daemi–Scaduto's S-complex
of singular instanton homology for the torus knots T(3,n), with the matching facts for two-bridge
knots. Each result is attributed to its source in the theorem headers: Lewallen,
Hedden–Herald–Kirk and Nagasato–Yamaguchi for the two-bridge characters (all binary dihedral),
Klassen for the dihedral count, and Daemi–Scaduto and Li–Ye for the gradings and ranks.

The paper records the consequence that, for every torus knot, the S-complex differential has rank
(1 + |σ| − ‖Δ‖₁)/2. For T(3,n) this is 0 when n ≡ 1, 2 and 1 when n ≡ 4, 5 (mod 6), for all n.
It also gives a Khovanov-homology proof that rank I♮(T(3,n)) = ‖Δ‖₁, using Schütz's
decomposition of three-braid complexes, and it compares the rank-one differential of 8₁₉ = T(3,4)
with the pillowcase bigon of Hedden–Herald–Kirk. The numerical checks are `riley_check.py`,
`torus_characters.py` and `fs_gradings.py` in `pillowcase/`.

The current version cites papers 2, 4 and 5 where they continue this material. arXiv:2607.26095 has
an earlier version.

### Paper 4: the S-complex differential of torus knots (`paper4/`)

For a torus knot K = T(p,q), the S-complex built on the traceless character variety has rank
1 + |σ(K)|, and its homology I♮(K) has rank ‖Δ_K‖₁. So the differential has rank
R(p,q) = (1 + |σ| − ‖Δ‖₁)/2. The paper proves:

- **A floor-sum formula** for R(p,q).
- **A sharp lower bound,** R(p,q) ≥ ⌊(p − 2b)²/4⌋ for coprime 2 ≤ p < q, where b is the inverse
  of q mod p.
- **A perfectness classification.** The differential vanishes exactly when p = 2, or when p is odd
  and q ≡ 2 or q ≡ p + 2(−1)^((p−1)/2) (mod 2p). The "if" direction is due to Daemi and Scaduto.
- **The graded components.** The components from degree 3 to 2 and from degree 1 to 0 have ranks
  ⌊R/2⌋ and ⌈R/2⌉, independently of the choices in the construction. So the ℤ/4-graded group
  I♮(T(p,q); ℚ) depends only on σ and ‖Δ‖₁, and the paper writes it down explicitly.
- **A comparison with the pillowcase.** In each of Hedden–Herald–Kirk's computations, a bigon ends
  at the reducible generator exactly when Daemi–Scaduto's invariant h is positive.
- **A conjecture** identifying the two complexes at chain level over 𝔽₂; paper 5 proves it for T(3,n).

### Paper 5: pillowcase Floer homology of T(3,n) (`paper5/`)

Hedden, Herald and Kirk defined a Lagrangian Floer homology in the pillowcase from a tangle
decomposition of a knot, and conjectured that for a suitable decomposition it recovers I♮. Paper 5
computes this homology for **every** torus knot T(3,n), using their decomposition:

- **The complex.** For all small holonomy perturbations in an explicit open cone, the immersed curves
  are restricted and the complex has 1 + |σ| generators.
- **The differential.** There is none when n ≡ 1, 2 (mod 6), and exactly one, to the abelian
  generator, when n ≡ 4, 5 (mod 6).
- **The homology** has the ℤ/4-graded dimensions of I♮.

This verifies the Hedden–Herald–Kirk conjecture for an infinite family with nonzero differentials.
Combined with a cobordism inequality of Daemi and Scaduto, it also proves paper 4's chain-level
conjecture for T(3,n): over 𝔽₂ the pillowcase complex is S-chain homotopy equivalent to the instanton
S-complex.

The proof rests on an explicit description of the traceless character variety and on the strict
monotonicity of an angle function along its components. Three inequalities are verified by
interval arithmetic. Section 7 reports a numerical computation, which is not a proof: for a
different decomposition of T(3,4) in the same family, the homology with zero bounding cochain has
rank 7 rather than 5. The source is `paper5/main.tex` with `paper5/sections/*.tex`.

The computations are in [`paper5/code/`](paper5/code). The three inequalities of Appendix A are
proved there by interval arithmetic (`mpmath`), and the numerical computations of Sections 5–7 are
reproduced. Each script prints every number the paper states, next to the computed value, and the
whole set runs in about twelve minutes (`cd paper5/code && ./run_all.sh`; it needs `numpy` and
`mpmath`, see `requirements.txt`).

Paper 2 builds on this paper: shearing its curves gives the pillowcase homology of the pretzel knots
P(−2,3,q) and of a family of twisted torus knots.

### Paper 2: the (−2,3,q) pretzel knots (`paper2/`)

The pretzel knots P(−2,3,q), q odd, begin with the torus knots T(3,4) and T(3,5) and are hyperbolic
for q ≥ 7. The paper studies both sides of the correspondence for this family.

- **Instanton homology.** For every odd q ≥ 3, I♮(P(−2,3,q); ℤ) is free abelian of rank q + 2. Over
  ℚ this is due to Lobb and Zentner; for q ≥ 7 the integral and graded statements also follow from
  work of Daemi and Scaduto. The proof here is uniform in q.
- **The conjecture of Hedden, Herald and Kirk.** Regluing the tangle that Hedden, Herald and Kirk use
  for T(3,5) by (q − 5)/2 Dehn twists along the Conway sphere gives a decomposition of P(−2,3,q).
  For small holonomy perturbations its pillowcase complex is computed by hand: it has q + 2 + 4N_q
  generators, a differential of rank 2N_q (nonzero for q ≥ 9), and homology isomorphic to I♮ as a
  ℤ/4-graded vector space. This verifies the conjecture for every knot in the family. The same
  argument applies to the twisted torus knots T(3,n;2,m) with n ≡ 1, 2 (mod 6).
- **Smith's decomposition.** For Smith's decomposition into Q₋₁/₂ and Q₁/₃ + Q₁/q, the pairing with
  zero bounding cochain has the wrong rank. In the two-arc algebra of Kotelskiy, Watson and
  Zibrowius the Maurer–Cartan equation is a finite identity. At q = 7, smoothings of a model of
  Smith's curve give Maurer–Cartan elements with the correct pairing, and a second closure selects
  one homotopy class among them, conditionally on the instanton–pillowcase correspondence and a
  stated localization hypothesis. The Lagrangian correspondence induced by the line of Q₁/₃ is
  immersed with a triple point, so Gao's representability theorem does not apply to it.
- **The Khovanov arc.** For q = 5, 7 the selected object is the sum of a rational arc and a closed
  curve with a two-dimensional local system, and one further Maurer–Cartan term near a corner turns it
  into the Bar-Natan invariant of the tangle; no rational closure distinguishes the two. The paper
  asks, without conjecturing an answer, whether the instanton object of the tangle is its Bar-Natan
  invariant.

The source is `paper2/main.tex` with `paper2/sections/*.tex`. The computations are in
[`paper2/code/`](paper2/code) (Sections 4, 5 and 7, Appendix A) and [`pillowcase/`](pillowcase)
(Section 6); see the README files there.

## Building

Build any paper with `pdflatex main.tex` from its directory (run it twice for cross-references).
Each paper compiles on its own with a standard TeX distribution. Compiled PDFs are included under
distinctive names: `paper1/traceless-gradings.pdf`, `paper2/pretzel-instanton-pillowcase.pdf`,
`paper4/s-complex-torus-knots.pdf` and `paper5/pillowcase-torus-knots.pdf`.

## The code

| directory | papers | requirements |
|---|---|---|
| [`pillowcase/`](pillowcase) | paper 2, Section 6; the numerical checks of paper 1 | Python 3, standard library |
| [`paper2/code/`](paper2/code) | paper 2, Sections 4, 5 and 7 and Appendix A | Python ≥ 3.11; numpy, sympy, mpmath for Sections 4–5 |
| [`paper5/code/`](paper5/code) | paper 5 | Python ≥ 3.11; numpy, mpmath |

Every program prints each number its paper states, next to the computed value, with `PASS` or
`FAIL`, and exits with status 0 if and only if all its checks pass. Each directory has a README
listing what each program checks, and a script `run_all.sh`:

```bash
sh pillowcase/run_all.sh fast            # paper 2, Section 6, without the slow finite tables
cd paper2/code && ./run_all.sh --quick   # paper 2, Sections 4, 5 and 7 (about 3 minutes)
cd paper5/code && ./run_all.sh           # paper 5 (about 12 minutes)
```

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

Code is released under the MIT License (`LICENSE`). The papers (`paper1/`, `paper2/`, `paper4/`, `paper5/`) are
© Bernd J. Wuebben; you may read and redistribute them for scholarly purposes with attribution.
