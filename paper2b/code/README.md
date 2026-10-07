# Programs for Section 4, Remark 3.17 and Appendix A

These programs carry out the computations of Section 4, Remark 3.17 and Appendix A of
*Maurer–Cartan elements in the pillowcase and the Bar-Natan invariant for the (−2,3,q) pretzel knots*. Statement
numbers refer to the compiled paper. The programs for the rest of Section 3 are in the directory `pillowcase/`.

## Smith's decomposition and the Khovanov arc

These scripts reproduce the computations reported in Section 4, in Remark 3.17 and in Appendix A.1, A.2 and A.4, and recheck
the morphism dimensions of Propositions 3.7 and 3.8 with two independent implementations. Statement numbers refer to
the compiled paper. There are two kinds of computation.

* **Exact computations over F₂.** Twisted complexes over the algebra 𝓑 of Kotelskiy, Watson and Zibrowius are
  finite matrices, the Maurer–Cartan equation is the finite identity (δ + b)² = 0, strict isomorphism of loop-type
  complexes is decided by comparing the words of their components, and morphism homology is computed in full, without
  truncating word length. These computations are the proofs of Theorem 4.1 and of the table in the proof of
  Corollary 4.6, and they establish the statements of Remarks 4.2, 4.3, 4.7, 4.8, 3.17, Proposition 4.10 and
  Appendix A.1–A.2. Their inputs are the complexes printed in the paper, curves encoded by `kwz_encode.py`, the
  invariants of Kotelskiy, Watson and Zibrowius computed by Zibrowius's program kht++, and ranks of reduced Khovanov
  homology computed with Khoca (both shipped as data, see below). The searches of Remarks 4.2 and 4.11 are exhaustive
  within the classes they state. Theorem 4.5 is proved in the paper by a geometric argument; the census of linear
  curves in `kwz_linear_pairings.py` supports Remark 4.7.
* **Floating-point input.** The encodings of the numerical traces of Smith's perturbed curve (Appendix A.4, "Traced
  curves") come from floating-point computations, and the statements about them (Remark 4.11, first sentences, and
  the cases q = 17, 19 of its searches) are corroboration only.

Every script prints each number the paper states next to the computed value, with `PASS` or `FAIL` (checks of
consistency that the paper does not state are marked `expected` instead of `paper`), and exits with status 0 if and
only if all its checks pass.

| script | verifies | runtime |
|---|---|---|
| `kwz_encoder_check.py` | Appendix A.4, "Encodings": for each of the 216 rational tangles Q_s computed by kht++, BN~(Q_s) is the straight lattice arc and Kh~(Q_s) the figure-eight around it. Section 4.1: Φ(Kh~(Q_{−r})) is the earring E_r at 49 slopes (and at all 342 slopes of the set of Corollary 4.6). The encodings of E_{−1/2}, E_{−3/4} are the complexes (3.2), (3.3) | 2 s |
| `kwz_hom_check.py` | Propositions 3.7 and 3.8 and Appendix A.4, "Morphism homology": the dimensions 7, 25; 9, 31; 9, 31; 9, 23; 9, 25 by both methods; the data (e, τ, r) = (228, 251, 227), (228, 251, 224), (228, 251, 228) printed after (3.3); f_L = f_M = f_R = 0 and p_{X,i} = U_X; Remark 3.6 (End(N) has f_M = f_R = 1, End(ι_∘, 0) has f_L = f_M = 1); agreement of the two methods on 34 pairings of Section 4 | < 1 s |
| `bn_structure.py` | Theorem 4.1 (a)–(c): (N_q, b) ≅ α_q ⊕ (c, J₂) with components of 23 + 8 and 15 + 8 generators; BN_q is the slide of α_q at (0,0) through two periods, the only one of the 24 slides, for q = 5, 7, …, 21; the Maurer–Cartan identities for κ_q and (N_q, b + κ_q) ≅ BN_q. Remark 4.3: the slide at (π, 0), (N_q, b + κ′_q) ≅ BN~(T_q), and the pairings 13 against 11 (q = 7) and 9 against 7 (q = 5) with E_{−1}. Appendix A.1: N_5, its ends 0 and 22, its pairings 5 and 15, b_5, the components of (N_5, b_5), the arc of (N_5, b_5 + κ_5) and its pairings 7 and 21. Also N_7 and N_5 against the encodings of the piecewise-linear models. Section 4.2: c and (c, J₂) admit a δ-grading but no bigrading, α_q and BN_q are arcs | < 1 s |
| `corner_terms.py` | Remark 4.2 and Appendix A.2: the one-arrow replacements of (N_q, b) (566 candidates for each b at q = 7, 402 at q = 5), exactly two giving BN_q in each case, as tabulated, with the added arrow labelled R at the end over (0,0); κ₇ common to b_{S18} and b_{S25}; the terms for t > 0. Edit distance 2 from N_q to BN_q with two minimal edits, at q = 7 both R-arrows touching 36. The 110 and 58 algebraic switches and the 221 925 and 32 567 sums of at most three, none ≅ BN_q. The census switches at q = 7 against the algebraic switches | 16 s |
| `closures_rational.py` | (4.1) at the 340 (q = 5) and 342 (q = 7) slopes of Corollary 4.6; the three numbers of Corollary 4.6 agree at every slope; its table at E_{−1} and E_{r_q} (7, 30; 11, 46) and the undeformed values 7, 28 and 11, 44; 108 slopes with odd numerator and denominator; N_q below rk Kh~ at 40 and 43 slopes and never above; the values q + 2 at r = −1/2 and 31 at r = −3/4. `--khoca` recomputes the Khovanov ranks with Khoca | 19 s |
| `kwz_linear_pairings.py` | Remark 4.7: 1887 linear curves for each q = 5, 7, …, 13, of which exactly 22 pair differently with BN_q and with α_q ⊕ (c, J₂), all parallel to c (the 19 side patterns on the line through (0,0) and (π,π), the two constant patterns on the line through (π,0) and (0,π), the simple closed curve), none parallel to α_q; the values for c and for the figure-eight around the arc from (0,0) to (π,π). Remark 4.8: 11/13, 22/26, 4, 52/44 (q = 7) and 7/9, 14/18, 36/28 (q = 5). Proposition 4.10. `--quick`: Remark 4.7 for q = 5, 7 only | 2.5 min (`--quick` 49 s) |
| `sign_positive_t.py` | Remark 3.17: for t > 0 the undeformed pairings 7, 25; (9, 31) at S₁₈, S₂₅, (5, 25) at S₆₉ and (7, 25) at S₇₄; the seven switches reaching nine; among the 82 census smoothings, the 82 opposite smoothings and the 16 span elements, (9, 31) only at S₁₈ and S₂₅; the third closure 69 against 103; in the span only b_{S18}, b_{S25} reach nine. The same quantities for t < 0 | 1 s |
| `yoneda_example.py` | Example 3.19: (N, b_{S18}) and (N, b_{S18} + b_{S35}) are Maurer–Cartan with the pair (9, 31); the figure-eight around the line of slope 1/3 through (0,0) and (π,π) pairs to 9 and 11 with them; the one around the parallel line through (π,0) and (0,π) pairs to 9 with both | < 1 s |
| `kwz_traced_curves.py` | Appendix A.4, "Traced curves", and Remark 4.11, first sentences (floating-point input): for q ≤ 13 the traces at a given q have strictly isomorphic encodings, isomorphic to N_q, with 47 and 55 generators at q = 11, 13; the pairing with E_{−1/2} is q (q = 5, 7), q + 4 (q = 11, 13), q + 8 (q = 17, 19); one switch raises it to q + 2 at q = 5, 7 | < 1 s |
| `switches_search.py` | Remark 4.11: edit distances 2, 2, 4, 4 and at least 6 (q = 17, 19, edit distances ≤ 5 excluded); chains of algebraic switches followed by one replacement: one switch at q = 5, 7, none with at most three at q = 11, 13, the first with four at q = 11, ending at an arc plus a closed curve other than α₁₁ ⊕ (c, J₂); switch distance from N_q to α_q ⊕ (c, J₂) 1 at q = 5, 7 and more than 5 at q = 11; the 35 004 complexes (6 331 up to strict isomorphism) obtained from N₁₁ by at most two switches: none has α₁₁ as a summand and none is compatible with the instanton data at the first closure and on the panel of 25 slopes. `--quick`: without the four-switch search, the switch distance at q = 11, the cases q = 17, 19, and with 2000 of the complexes | 30 min (`--quick` 40 s) |

Shared code:

* `kwz_algebra.py`: the algebra, twisted complexes in the notation of the paper, the parser for the output of kht++,
  words and strict isomorphism of loop-type complexes, and the two implementations of the morphism homology: the
  mapping-cone reduction of Lemma 3.5 (the free ranks f_X and the torsion are read off from quotients by powers of
  U_X, and the rank of β_* is computed in such a quotient) and the structure of Mor(X, Y) as a free complex over
  F₂[H], H = D + S² central.
* `kwz_encode.py`: the encoding of curves in the pillowcase as twisted complexes (planar cover, arc system, Smith's
  coordinates with the corner (0, π) at infinity, the translation Φ), straight arcs, figure-eights, linear curves, the
  earrings E_r, KWZ's rational invariants, α_q, c, the doubled word of c, BN_q and the 24 slides.
* `kwz_objects.py`: the complexes printed in the paper (N_7 = N, E, E_*, the switches of Section 3.3, N_5, b_5, κ_q,
  κ′_q), the slope set of Corollary 4.6, and loaders for the data files.
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

Python 3.11 or later; the programs use only the standard library. Khoca and SnapPy are needed only for
`closures_rational.py --khoca`. Tested with Python 3.11.14 and 3.14.4.

```bash
./run_all.sh            # everything, about 35 minutes (PY=/path/to/python to choose the interpreter)
./run_all.sh --quick    # shorter runs of kwz_linear_pairings.py and switches_search.py, about 2 minutes
```
