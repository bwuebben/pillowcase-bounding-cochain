#!/usr/bin/env python3
r"""
check_paper1_numbers.py -- numeric-integrity gate for paper 1
("Traceless SU(2) characters and Z/4 instanton gradings for two-bridge and
(3,n)-torus knots").

Recomputes every published quantity from first principles (via the ancillary
scripts' own logic) and asserts each appears in the compiled PDF, anchored to
nearby text so that a wrong number cannot pass by matching a coincidental
occurrence elsewhere.  Then sweeps the remaining numerals of the body (the
bibliography is excluded) and lists any that no registered fact explains.

Run from pillowcase/ or paper/anc/:   python3 check_paper1_numbers.py
Requires pdftotext (poppler) and ../paper/main.pdf (or ./main.pdf).
Exit 1 on any registered-fact miss; the residue list is informational.
"""

import re
import subprocess
import sys
from math import pi, cos, sin, gcd, floor, sqrt, acos, degrees
from pathlib import Path

HERE = Path(__file__).resolve().parent
for cand in [HERE.parent / "paper" / "main.pdf", HERE / "main.pdf",
             HERE.parent / "main.pdf"]:
    if cand.exists():
        PDF = cand
        break
else:
    sys.exit("main.pdf not found")

sys.path.insert(0, str(HERE))
from riley_check import phi_exact                       # noqa: E402
from fs_gradings import valid, mu, gr_FS, rho_FS        # noqa: E402


def pdf_text():
    out = subprocess.run(["pdftotext", str(PDF), "-"],
                         capture_output=True, text=True, check=True)
    t = out.stdout
    body = t.split("\nReferences\n")[0]                 # exclude bibliography
    return re.sub(r"\s+", " ", body)                    # normalized whitespace


BODY = pdf_text()
failures = []
registered_numbers = set()


def check(desc, needle, anchor=None, window=300):
    """Assert needle appears in the body; if anchor given, within `window`
    chars of it (presence-in-context for non-unique values)."""
    needle_n = re.sub(r"\s+", " ", needle)
    ok = False
    if anchor is None:
        ok = needle_n in BODY
    else:
        a = re.sub(r"\s+", " ", anchor)
        start = 0
        while True:
            i = BODY.find(a, start)
            if i < 0:
                break
            lo, hi = max(0, i - window), i + len(a) + window
            if needle_n in BODY[lo:hi]:
                ok = True
                break
            start = i + 1
    print(("  PASS " if ok else "  FAIL ") + desc)
    if not ok:
        failures.append(desc)
    for tok in re.findall(r"\d+", needle_n):
        registered_numbers.add(tok)


def fmt_poly(coeffs):
    """Integer poly (index=degree) -> pdftotext-style string, descending."""
    terms = []
    for d in range(len(coeffs) - 1, -1, -1):
        c = coeffs[d]
        if c == 0:
            continue
        mag = "" if (abs(c) == 1 and d > 0) else str(abs(c))
        var = "" if d == 0 else ("u" if d == 1 else f"u{d}")
        terms.append(("+ " if c > 0 else "- ") + mag + var)
    s = " ".join(terms)
    return s[2:] if s.startswith("+ ") else s


# ---------------------------------------------------------------------------
# recomputation layer
# ---------------------------------------------------------------------------

def N_formula(n):
    r = n % 6
    return {1: 2 * (n - 1) // 3, 2: (2 * n - 1) // 3,
            4: (2 * n + 1) // 3, 5: 2 * (n + 1) // 3}[r]


def meridian(n):
    for a in range(-3, 4):
        if (1 - a * n) % 3 == 0:
            return a, (1 - a * n) // 3


def N_enum(n):
    a, b = meridian(n)
    cnt = 0
    for l1 in (1, 2):
        for l2 in range(1, n):
            if (l1 - l2) % 2:
                continue
            c = (cos(a * pi * l1 / 3) / sin(a * pi * l1 / 3)) * \
                (cos(b * pi * l2 / n) / sin(b * pi * l2 / n))
            if -1 < c < 1:
                cnt += 1
    return cnt


def alex_sum(n):
    """Sum of |coefficients| of Delta_{T(3,n)} by exact polynomial division."""
    def cyc(N):
        return {N: 1, 0: -1}

    def mul(p, q):
        r = {}
        for i, ci in p.items():
            for j, cj in q.items():
                r[i + j] = r.get(i + j, 0) + ci * cj
        return {k: v for k, v in r.items() if v}

    def div(num, den):
        num = dict(num)
        q = {}
        dd, dc = max(den), den[max(den)]
        while num:
            nd = max(num)
            if nd < dd:
                assert not any(num.values())
                break
            c = num[nd] // dc
            q[nd - dd] = c
            for k, v in den.items():
                num[nd - dd + k] = num.get(nd - dd + k, 0) - c * v
                if num[nd - dd + k] == 0:
                    del num[nd - dd + k]
        return q
    d = div(div(mul(cyc(3 * n), cyc(1)), cyc(3)), cyc(n))
    return sum(abs(v) for v in d.values())


# ---------------------------------------------------------------------------
print(f"checking {PDF} against recomputed values\n")

# (1) the four printed Riley polynomials (Remark 3.1), from the exact recursion
for p in (3, 5, 7, 9):
    poly = fmt_poly(phi_exact(p))
    check(f"Remark 3.1: phi_{p} = {poly}", poly, anchor="Small cases")
    # constant term p and degree (p-1)/2 are implied by the string itself
    assert phi_exact(p)[0] == p and len(phi_exact(p)) - 1 == (p - 1) // 2

# (2) Prop 4.1 printed values, from formula AND enumeration
vals = [N_enum(n) for n in (2, 4, 5, 7, 8, 10, 11, 13)]
assert vals == [N_formula(n) for n in (2, 4, 5, 7, 8, 10, 11, 13)]
assert all(N_formula(n) == floor(5 * n / 6) - floor(n / 6) == N_enum(n)
           for n in range(2, 200) if gcd(3, n) == 1)   # the closed-form proof
check("Prop 4.1 values: " + ", ".join(map(str, vals)),
      ", ".join(map(str, vals)), anchor="are")

# (3) Example 5.1: Anvari calibration, raw values recomputed
for l3, want_gr, want_mod8, want_mu, want_mod4 in [(2, 175, 7, 87, 3),
                                                   (4, 251, 3, 125, 1)]:
    m, g, r = mu(7, 2, l3)
    raw = round(0.5 * g + 0.25 * (1 - r))
    assert (g, g % 8, raw, m) == (want_gr, want_mod8, want_mu, want_mod4), \
        (l3, g, raw, m)
    check(f"Ex 5.1: gr = {want_gr} = {want_mod8} (mod 8)",
          f"gr = {want_gr} ≡ {want_mod8}", anchor="(1, 2, ")
    check(f"Ex 5.1: mu = {want_mu} = {want_mod4} (mod 4)",
          f"µ = {want_mu} ≡ {want_mod4}", anchor="(mod 8)")

# (4) Table 1, all four numeric rows recomputed
ns = (5, 7, 11, 13, 17, 19, 23, 25)
N_row = [N_formula(n) for n in ns]
a_row = [N // 2 for N in N_row]
chain_row = [1 + 4 * a for a in a_row]
inat_row = [alex_sum(n) for n in ns]
assert inat_row == [(1 + 4 * a) if n % 6 == 1 else (4 * a - 1)
                    for n, a in zip(ns, a_row)]
# triangle-count agreement (eq (2)) for the full claimed range n <= 43
for n in range(5, 44, 2):
    if gcd(3, n) == 1:
        A = len([l3 for l3 in range(2, n, 2) if valid(n, l3)])
        assert N_formula(n) == 2 * A, n
        # and the even mu split in {1,3} (Prop 5.2)
        mus = [mu(n, 2, l3)[0] for l3 in range(2, n, 2) if valid(n, l3)]
        assert sorted(set(mus)) in ([1, 3], [1], [3]) and \
            mus.count(1) == mus.count(3), n
# pdftotext linearizes the table in column blocks, so rows cannot be matched
# as strings; instead compare the MULTISET of numerals in the table window
# (between the last header cell and the caption) against all 48 recomputed
# entries.  A single wrong cell changes the multiset and fails the check.
mods_row = [n % 6 for n in ns]
expected = sorted(map(str, list(ns) + mods_row + N_row + a_row
                      + chain_row + inat_row))
hdr = BODY.find("chain rank 1 + 4a rank I ♮")
cap = BODY.find("Table 1.", hdr)
if hdr < 0 or cap < 0:
    print("  FAIL Table 1: window not found")
    failures.append("Table 1 window")
else:
    got = sorted(re.findall(r"\d+", BODY[hdr + len("chain rank 1 + 4a rank I ♮"):cap]))
    ok = got == expected
    print(("  PASS " if ok else "  FAIL ")
          + f"Table 1: all 48 cells match recomputation (multiset)")
    if not ok:
        failures.append("Table 1 multiset")
        print(f"        expected {expected}\n        got      {got}")
    registered_numbers.update(expected)

# (5) 8_19: the three theta values, complexes, ranks
a, b = meridian(4)
th = {}
for l1, l2 in [(1, 1), (1, 3), (2, 2)]:
    c = (cos(a * pi * l1 / 3) / sin(a * pi * l1 / 3)) * \
        (cos(b * pi * l2 / 4) / sin(b * pi * l2 / 4))
    th[(l1, l2)] = degrees(acos(c))
assert abs(th[(1, 1)] - 125.26) < 0.05 and abs(th[(1, 3)] - 54.74) < 0.05 \
    and abs(th[(2, 2)] - 90) < 1e-9
check("8_19 thetas 125.3, 54.7, 90", "125.3◦ , 54.7◦", anchor="non-dihedral")
check("8_19 dihedral at 90", "90◦", anchor="one dihedral")
check("8_19 chain (2,1,2,2)", "(2, 1, 2, 2)", anchor="χ = +1")
check("8_19 homology (2,1,1,1)", "(2, 1, 1, 1)", anchor="I ♮ (819")
check("8_19 rank 5 > det 3", "5 > det = 3", anchor="non-thin")

# (6) headline abstract facts and knot names
assert alex_sum(5) == 7
check("abstract: T(3,5) rank 7 not 9", "rank 7, not 9", anchor="10124")
check("knot name 10_124 for T(3,5)=P(-2,3,5)", "10124", anchor="P (−2, 3, 5)")
check("knot name 8_19 = T(3,4)", "819 = T (3, 4)")

# (7) verification-range claims match what the shipped scripts do
check("claim: enumeration for all n <= 25", "n ≤ 25", anchor="direct enumeration")
check("claim: eq (2) verified odd n <= 43", "n ≤ 43",
      anchor="triangle-inequality count")
check("claim: gradings evaluated odd n <= 43", "n ≤ 43",
      anchor="direct evaluation")
check("claim: Riley words for all p <= 9", "p ≤ 9",
      anchor="Gaussian-integer")

# (8) residue sweep: numerals in the body no registered fact explains
STRUCTURAL = set("0 1 2 3 4 5 6 7 8 9 10 11 12 13".split())  # small ints:
# section/case labels, small indices, residues mod 6, footnote-scale ints
years = {str(y) for y in range(1984, 2027)}
msc = {"57R58", "57K10", "57K31", "53D40", "58J30", "2020"}
residue = []
for m in re.finditer(r"\d+(?:\.\d+)?", BODY):
    tok = m.group(0)
    if tok in registered_numbers or tok in STRUCTURAL or tok in years \
            or tok in msc:
        continue
    if re.fullmatch(r"\d\.\d\d?", tok):        # section / result labels
        continue
    pre = BODY[max(0, m.start() - 8):m.start()]
    if "[" in pre or re.search(r"[A-Z][a-z]*$", pre):   # citation labels
        continue
    ctx = BODY[max(0, m.start() - 30):m.end() + 30]
    residue.append((tok, ctx.strip()))
print(f"\nresidue: {len(residue)} unexplained numerals (informational)")
for tok, ctx in residue[:40]:
    print(f"    {tok:>8}  …{ctx}…")

print()
if failures:
    print(f"FAILURES: {len(failures)}")
    sys.exit(1)
print("ALL REGISTERED FACTS PRESENT IN PDF")
