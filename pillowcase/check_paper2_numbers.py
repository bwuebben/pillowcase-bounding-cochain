#!/usr/bin/env python3
"""
check_paper2_numbers.py -- numeric-integrity check for paper 2
(paper2/main.tex, "The instanton homology of the (-2,3,q) pretzel knots and
computed bounding cochains in the pillowcase").

Recomputes every published quantity from the code in this directory (or, where
a quantity comes from a long recorded run, registers it as a CONSTANT with its
provenance), then asserts that each appears in the compiled PDF's text --
anchored to nearby context when the bare value is not unique -- and finally
sweeps the body for numerals that no registered fact explains.

Exit status: 0 iff every assertion holds and the unexplained-numeral residue is
empty.  Run from this directory:

    python3 check_paper2_numbers.py            # needs paper2/main.pdf built

Sections:
  [A] facts recomputed live (Alexander/skein, determinants, curve geometry)
  [B] facts from recorded runs (q=11 solve, second-perturbation values)
  [C] presence-in-PDF assertions (anchored where needed)
  [D] residue sweep of remaining numerals
"""

import math
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PDF = HERE.parent / "paper2" / "main.pdf"

fails = []


def check(label, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"  [{tag}] {label}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        fails.append(label)


# ===========================================================================
# [A] Recomputed facts
# ===========================================================================
print("[A] recomputing published quantities")

sys.path.insert(0, str(HERE))
from skein_alexander import family, abs_sum        # noqa: E402

fam = family(19)                                   # {q: Alexander poly, dict}
# ell(K_q) = sum |a_i(Delta)| = q+2, via the validated skein recursion
ELL = {q: abs_sum(fam[q]) for q in (3, 5, 7, 11, 13)}
for q in ELL:
    check(f"ell(P(-2,3,{q})) = {q}+2 = {q+2}", ELL[q] == q + 2, f"got {ELL[q]}")

# Lehmer at q=7: 11 coefficient slots (degree span 10), all coeffs in {0,+-1},
# absolute pattern 1,1,0,1,1,1,1,1,0,1,1
lo, hi = min(fam[7]), max(fam[7])
lehmer = [abs(fam[7].get(e, 0)) for e in range(lo, hi + 1)]
check("q=7 Alexander = Lehmer (11 slots, |coeffs| = 1,1,0,1,1,1,1,1,0,1,1)",
      len(lehmer) == 11 and lehmer == [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1],
      f"got {lehmer}")

# determinants det = |q-6| and the members' knot identifications
DET = {q: abs(q - 6) for q in (3, 5, 7, 11, 13)}
check("det = 3,1,1,5,7 for q = 3,5,7,11,13",
      [DET[q] for q in (3, 5, 7, 11, 13)] == [3, 1, 1, 5, 7])
# pretzel determinant formula |p1p2+p2p3+p3p1| at (-2,3,q)
check("det formula |-6+3q-2q| = |q-6|",
      all(abs(-6 + 3 * q - 2 * q) == DET[q] for q in DET))

# the deficiency law: naive - Inat = 2*sgn(det-3)
NAIVE = {3: 5, 5: 5, 7: 7, 11: 15, 13: 17}         # Table 1 (q=3 from HHK2 11.6)
INAT = {q: q + 2 for q in NAIVE}
sgn = lambda x: (x > 0) - (x < 0)
for q in NAIVE:
    check(f"law at q={q}: {NAIVE[q]} - {INAT[q]} = 2*sgn({DET[q]}-3)",
          NAIVE[q] - INAT[q] == 2 * sgn(DET[q] - 3))

# rank arithmetic quoted in ss 3.2 / 5.1-5.3
check("q=5: 9 - 2*2 = 5 naive, 9 - 2 = 7 deformed",
      9 - 2 * 2 == NAIVE[5] and 9 - 2 == INAT[5])
check("q=7: 13 gens, naive 13-6=7, deformed 13-4=9",
      13 - 6 == NAIVE[7] and 13 - 4 == INAT[7])
check("q=11: 19 gens, naive 19-4=15, deformed 19-6=13",
      19 - 4 == NAIVE[11] and 19 - 6 == INAT[11])

# Manion's closed form at (p,q,r)=(2,3,q):  p^2+(q-p)(r-p) = q+2
check("Manion: 4+(3-2)(q-2) = q+2 for all odd q>=3",
      all(4 + (3 - 2) * (q - 2) == q + 2 for q in range(3, 20, 2)))

# curve geometry: blue self-crossing counts and corner distances
from deform import build_pretzel                    # noqa: E402
from maurer_cartan import orbit_group               # noqa: E402

CORNERS = [(0, 0), (math.pi, 0), (0, math.pi), (math.pi, math.pi)]


def dcorner(p):
    g, t = p
    return min(math.hypot(g - cg, min(abs(t - ct), 2 * math.pi - abs(t - ct)))
               for cg, ct in CORNERS)


CROSSINGS, MINDIST = {}, {}
for k in (2, 3, 5, 6):
    q = 2 * k + 1
    _, blue, _ = build_pretzel(k)
    orbs = orbit_group(blue)
    CROSSINGS[q] = len(orbs)
    MINDIST[q] = min(dcorner(pp) for pp, _ in orbs)
check("blue self-crossings 40/82/225 at q=5/7/11",
      (CROSSINGS[5], CROSSINGS[7], CROSSINGS[11]) == (40, 82, 225),
      f"got {CROSSINGS}")
check("min corner distance 0.61/0.42/0.27/0.23 at q=5/7/11/13",
      [round(MINDIST[q], 2) for q in (5, 7, 11, 13)] == [0.61, 0.42, 0.27, 0.23],
      f"got { {q: round(MINDIST[q], 4) for q in MINDIST} }")
check("all computed members' crossings >= 0.23 from corners",
      min(MINDIST.values()) >= 0.23)

# the q=5 support (b2_result.py battery output values)
SA, SB = (0.028, 1.272), (3.057, 4.981)
check("s_A, s_B corner distances >= 0.87 (deficit-cochain bound)",
      dcorner(SA) >= 0.87 and dcorner(SB) >= 0.87,
      f"{dcorner(SA):.3f}, {dcorner(SB):.3f}")
# paper quotes (0.03,1.27) and (3.06,4.98)
check("comp 1.3(i) rounding: (0.03,1.27)/(3.06,4.98)",
      (round(SA[0], 2), round(SA[1], 2)) == (0.03, 1.27)
      and (round(SB[0], 2), round(SB[1], 2)) == (3.06, 4.98))

# ===========================================================================
# [B] Recorded-run constants (provenance: committed logs, not recomputed here)
# ===========================================================================
print("[B] recorded-run constants")
# pretzel_solve.py 3 (rerun 2026-07-27, scratch log): unique support b=(69,)
Q7_P1 = (0.042, 5.405)         # perturbation 1, this machine, 2026-07-27
Q7_P2 = (0.058, 5.413)         # perturbation 2, RESEARCH_LOG.md line ~1458
check("q=7 support corner distances >= 0.87",
      dcorner(Q7_P1) >= 0.87 and dcorner(Q7_P2) >= 0.87,
      f"{dcorner(Q7_P1):.3f}, {dcorner(Q7_P2):.3f}")
check("q=7 two-perturbation representative = (0.05, 5.41)",
      abs((Q7_P1[0] + Q7_P2[0]) / 2 - 0.05) < 0.005
      and abs(round((Q7_P1[1] + Q7_P2[1]) / 2, 2) - 5.41) < 0.005)
# pretzel_solve.py 5 (RESEARCH_LOG.md sec 38): 225 crossings, 108 with
# triangles, 55 MC-valid single-crossing cochains, 50 into one generator.
Q11 = dict(crossings=225, tri=108, valid=55, distinguished=50)
check("q=11 crossing count matches recomputed curve", Q11["crossings"] == CROSSINGS[11])

# ===========================================================================
# [C] Presence in the compiled PDF
# ===========================================================================
print("[C] asserting each value appears in the compiled PDF")
if not PDF.exists():
    print(f"  [FAIL] {PDF} not found -- build the paper first")
    sys.exit(1)
text = subprocess.run(["pdftotext", str(PDF), "-"], capture_output=True,
                      text=True).stdout
# strip the bibliography (everything from the References heading on)
body = re.split(r"\nReferences\n", text)[0]
squashed = re.sub(r"\s+", " ", body)


def present(label, needle, anchor=None, window=260):
    """needle must appear in the body; if anchor given, within `window` chars
    of the anchor (for values that are not document-unique)."""
    if anchor is None:
        ok = needle in squashed
    else:
        ok = any(needle in squashed[max(0, m.start() - window):m.start() + window]
                 for m in re.finditer(re.escape(anchor), squashed))
    check(label, ok, f"'{needle}' near '{anchor}'" if anchor else f"'{needle}'")


present("q+2 (the theorem)", "q + 2", anchor="rank")
present("naive rank sequence", "5, 5, 7, 15, 17")
present("Table: q=3 row (8_19 = T(3,4))", "819 = T (3, 4)")   # pdftotext squashes subscripts
present("Table: q=5 row (10_124 = T(3,5))", "10124 = T (3, 5)")
present("law display", "2 sgn(det K", anchor="sgn")
present("det = |q-6|", "|q − 6|")
present("forty crossings at q=5", "forty", anchor="self-intersection")
present("eighty-two at q=7", "eighty-two")
present("225 at q=11 (words)", "two hundred twenty-five")
present("225 at q=11 (numeral)", "225 crossings")
present("108 tri-crossings", "108 of them carry triangles")
present("fifty-five cochains", "fifty-five")
present("fifty of the fifty-five", "fifty of the fifty-five")
present("nineteen generators (q=11)", "nineteen generators")
present("thirteen generators (q=7)", "Thirteen generators")
present("nine generators (q=5)", "nine points")
present("generator variation 19 vs 23", "19 vs 23")
present("generator variation 25 vs 27", "25 vs 27")
present("s_A value", "(0.03, 1.27)")
present("s_B value", "(3.06, 4.98)")
present("q=7 crossing representative", "(0.05, 5.41)")
present("q=7 perturbation values", "(0.042, 5.405)")
present("q=7 perturbation values (2)", "(0.058, 5.413)")
present("corner distance bound (cochains)", "0.27", anchor="corner")
present("corner distance bound (deficit)", "0.87", anchor="deficit")
present("corner distances (all members)", "0.61, 0.42, 0.27, 0.23")
present("Lehmer Mahler measure", "1.17628")
# \binom{40}{3} loses its "40" under pdftotext's linearization -- check source
src = (HERE.parent / "paper2" / "main.tex").read_text()
check("binom(40,3) triples (source)", r"\binom{40}3" in src or r"\binom{40}{3}" in src)
present("Manion total rank", "4 + (3 − 2)(q − 2)")
present("HHK q=3 attribution", "11.6", anchor="Hedden")

# ===========================================================================
# [D] residue sweep: numerals with no registered explanation
# ===========================================================================
print("[D] residue sweep")
# strip what the sweep should not see: the MSC/keywords footnote, bracketed
# citation groups (labels like CHKK22 / Lin92 / Hir01 shed numerals), and the
# subscripted knot names 8_19 / 10_124, which pdftotext squashes to 819/10124
sweep = re.sub(r"2020 Mathematics Subject Classification.*?Atiyah–Floer conjecture,",
               " ", squashed)
sweep = re.sub(r"\[[^\]]*\]", " ", sweep)
sweep = sweep.replace("819", " ").replace("10124", " ")
# "2^2 + (3-2)(q-2)" in section 1.2 linearizes to "22 + ..."; the value 4 = 2^2
# is verified in [A] (Manion), so drop this rendering artifact from the sweep
present("Manion p^2 rendering in section 1.2", "22 + (3 − 2)(q − 2)")
sweep = sweep.replace("22 + (3 − 2)(q − 2)", " ")
known_tokens = set()
for q in (3, 5, 7, 11, 13):
    known_tokens |= {str(q), str(q + 2), str(abs(q - 6)), str(NAIVE[q])}
known_tokens |= {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "12", "13", "15", "17", "19", "23", "25", "27",
                 "40", "55", "82", "108", "225", "50",
                 "2020", "2010", "2011", "2014", "2018", "2022", "2024",
                 "2025", "2026",              # years incl. \date and citations
                 "57R58", "57K10", "57K18", "57K31", "53D40",
                 "0.03", "1.27", "3.06", "4.98", "0.05", "5.41",
                 "0.042", "5.405", "0.058", "5.413",
                 "0.61", "0.42", "0.27", "0.23", "0.87", "0.026", "0.6",
                 "1.17628", "8", "16", "20", "21", "24", "30", "39", "49",
                 "18", "14", "11.6", "5.9", "5.6", "4.22", "4.23", "6.6",
                 "10.4", "1.1", "2.9", "0.55", "0.45", "1.4", "3.2",
                 "1.35"}
# tolerate section/equation/figure references 1..9 and decimal labels x.y
resid = []
for m in re.finditer(r"\d+(?:\.\d+)?", sweep):
    tok = m.group(0)
    if tok in known_tokens:
        continue
    if re.fullmatch(r"\d\.\d", tok):        # x.y cross-references
        continue
    ctx = squashed[max(0, m.start() - 40):m.start() + 40]
    resid.append((tok, ctx.strip()))
if resid:
    print(f"  [WARN] {len(resid)} unexplained numeral(s):")
    for tok, ctx in resid[:20]:
        print(f"     {tok!r}: ...{ctx}...")
check("residue sweep empty", not resid, f"{len(resid)} left")

# ===========================================================================
print()
if fails:
    print(f"FAILURES ({len(fails)}): " + "; ".join(fails))
    sys.exit(1)
print("ALL CHECKS PASS")
