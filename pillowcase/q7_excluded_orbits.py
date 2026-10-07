#!/usr/bin/env python3
"""Single smoothings of the q=7 model that change the cancelled module.

At each of the 82 self-intersection orbits of the piecewise-linear model
L_{7,PL} there are two iota-invariant smoothings: the census smoothing
(``pairing=1`` in ``q7_quilt_census.physical_switch_only``) and the opposite
smoothing (``pairing=2``).  ``q7_quilt_census.py --all-node-census`` pairs the
73 census smoothings that preserve the 31-generator module of N.  This program
treats the remaining cases.  Every smoothed curve is encoded as a standalone
twisted complex over the algebra B, one representative for each pair of
involution-related torus components, and paired with rational earrings by the
exact mapping-cone reduction ``q7_kwz.exact_red_blue_homology``, which raises an
error unless every face complex has torsion homology.

It checks:

  * the census smoothing is not generator-preserving exactly at the nine orbits
    S7, S8, S23, S24, S49, S52, S53, S61, S73; each gives a 31-generator
    complex, with dimensions against the earrings of slopes -1/2 and -3/4
    equal to (3,23) at S7, (5,23) at S8, S53, S73, and (7,25) at the other five;
  * among all 82 census smoothings the pair (9,31) occurs only at S18, S25;
  * among all 82 opposite smoothings the pair (9,31) occurs only at S18, S25;
    there the smoothed curve has one component in the pillowcase, its
    encoding has 21 generators, and its dimension against the earring of
    slope 3/4 is 69, while that of (N, b_S18) is 103;
  * the numerator closure of Q_{3/4}+Q_{1/3}+Q_{1/7} has determinant 103.

Run from the repository root or from ``pillowcase``:

    python3 pillowcase/q7_excluded_orbits.py
"""
from __future__ import annotations

import contextlib
import io
import os
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import q7_kwz as q7                      # noqa: E402
import q7_closure_probe as closure       # noqa: E402
import q7_quilt_census as census         # noqa: E402
from bigons import edges_of              # noqa: E402
from maurer_cartan import orbit_group    # noqa: E402
from surgery_check import iota_pairing, smooth_curve  # noqa: E402

EPS = 0.006
EXCLUDED = {
    "S7": (3, 23),
    "S8": (5, 23), "S53": (5, 23), "S73": (5, 23),
    "S23": (7, 25), "S24": (7, 25), "S49": (7, 25), "S52": (7, 25),
    "S61": (7, 25),
}
SELECTED = ["S18", "S25"]

FAILURES = []


def report(label, ok):
    print(("[PASS] " if ok else "[FAIL] ") + label)
    if not ok:
        FAILURES.append(label)


def quiet(function, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return function(*args, **kwargs)


def standalone(components):
    """Twisted complex of a smoothed torus curve, one copy per iota-orbit.

    An iota-invariant component projects to an arc of the pillowcase and is
    encoded by its geometric precurve; two components exchanged by iota
    project to one closed curve and one of them is encoded cyclically.
    Returns (data, pre, delta, number of pillowcase components) in the format
    expected by ``q7_kwz.exact_red_blue_homology``.
    """
    encoded = []
    groups = defaultdict(list)
    for index, component in enumerate(components):
        cdata = q7.encode_type_d(component)
        invariant = len(cdata["hits"]) != len(cdata["unique"])
        pre = (q7.geometric_precurve(cdata) if invariant
               else q7.cyclic_precurve(cdata))
        if pre["residue"]:
            raise AssertionError("a component differential does not square to zero")
        encoded.append((cdata, pre))
        groups[(invariant, census._component_signature(cdata))].append(index)
    selected = []
    for (invariant, _), indices in groups.items():
        if invariant:
            selected.extend(indices)
        elif len(indices) == 2:
            selected.append(indices[0])
        else:
            raise AssertionError("non-invariant components are not an iota pair")
    unique, vertices, delta = {}, [], []
    for index in selected:
        cdata, pre = encoded[index]
        for v in pre["vertices"]:
            unique[(index, v)] = {"arc": cdata["unique"][v]["arc"]}
            vertices.append((index, v))
        delta.extend(((index, a), tuple(w), (index, b)) for a, w, b in pre["delta"])
    if any(w[0] not in ("L", "M", "R") or w[1] < 1 for _, w, _ in delta):
        raise AssertionError("encoded complex is not reduced")
    return {"unique": unique}, {"vertices": vertices}, q7.collect_f2(delta), len(selected)


def _deg(a):
    return a.bit_length() - 1


def _divmod(a, b):
    """Division with remainder in F2[U], polynomials encoded as integers."""
    q, db = 0, _deg(b)
    while a and _deg(a) >= db:
        shift = _deg(a) - db
        q ^= 1 << shift
        a ^= b << shift
    return q, a


def _mul(a, b):
    return q7._poly_mul_f2(a, b)


def _diagonalize(matrix, size, vectors):
    """Diagonalize a square matrix over F2[U] by row and column operations.

    ``matrix`` maps (row, col) to a polynomial; ``vectors`` are cycles whose
    coordinates undergo the same row operations.  Returns the nonzero
    diagonal entries and the transformed vectors.
    """
    m = [[matrix.get((i, j), 0) for j in range(size)] for i in range(size)]
    vec = [[v.get(i, 0) for i in range(size)] for v in vectors]
    diagonal = []
    t = 0
    while t < size:
        best = None
        for i in range(t, size):
            for j in range(t, size):
                if m[i][j] and (best is None or _deg(m[i][j]) < best[0]):
                    best = (_deg(m[i][j]), i, j)
        if best is None:
            break
        _, i0, j0 = best
        while True:
            m[t], m[i0] = m[i0], m[t]
            for v in vec:
                v[t], v[i0] = v[i0], v[t]
            for row in m:
                row[t], row[j0] = row[j0], row[t]
            pivot = m[t][t]
            dirty = False
            for i in range(t + 1, size):
                if m[i][t]:
                    q, r = _divmod(m[i][t], pivot)
                    if q:
                        m[i] = [x ^ _mul(q, y) for x, y in zip(m[i], m[t])]
                        for v in vec:
                            v[i] ^= _mul(q, v[t])
                    dirty |= bool(r)
            for j in range(t + 1, size):
                if m[t][j]:
                    q, r = _divmod(m[t][j], pivot)
                    if q:
                        for row in m:
                            row[j] ^= _mul(q, row[t])
                    dirty |= bool(r)
            if not dirty:
                break
            candidates = ([(i, t) for i in range(t + 1, size) if m[i][t]]
                          + [(t, j) for j in range(t + 1, size) if m[t][j]])
            i0, j0 = min(candidates, key=lambda c: _deg(m[c[0]][c[1]]))
        diagonal.append(m[t][t])
        t += 1
    return diagonal, vec


def mapping_cone(earring, blue):
    """Exact wrapped morphism homology via the mapping-cone reduction.

    Each face complex is a finite free F2[U]-complex.  Unit entries are
    cancelled by Gaussian elimination, the residual complex is diagonalized,
    and the identity images are followed through both steps.  Raises an error
    if some face complex has a free homology summand.  Returns the dimension,
    the numbers (e, tau, r) and the set of invariant factors.
    """
    red_data, red_pre = earring
    data, pre, delta = blue[:3]
    identities = q7._bounded_hom_basis(
        red_data, red_pre["vertices"], data, pre, 0)
    tau, factors, columns = 0, set(), [[] for _ in identities]
    for face in ("L", "M", "R"):
        labels, matrix = q7._face_tower_complex(
            face, red_data, red_pre["vertices"], red_pre["delta"],
            data, pre, delta)
        vectors = q7._identity_tail_vectors(
            face, labels, identities, red_pre["delta"], delta)
        active, residual, vectors, _ = q7._cancel_unit_tower_pairs(
            labels, dict(matrix), vectors)
        order = {old: new for new, old in enumerate(sorted(active))}
        small = {(order[i], order[j]): c for (i, j), c in residual.items()}
        moved = [{order[i]: c for i, c in v.items()} for v in vectors]
        diagonal, moved = _diagonalize(small, len(order), moved)
        if 2 * len(diagonal) != len(order):
            raise AssertionError(f"face {face}: free homology summand")
        for v in moved:
            if any(v[len(diagonal):]):
                raise AssertionError(f"face {face}: identity image is not a cycle")
        for k, s in enumerate(diagonal):
            factors.add(s)
            tau += _deg(s)
            for col, v in enumerate(moved):
                rem = _divmod(v[k], s)[1]
                columns[col].extend(
                    (face, k, b) for b in range(_deg(s)) if rem >> b & 1)
    positions = {}
    for col in columns:
        for key in col:
            positions.setdefault(key, len(positions))
    r = q7.sparse_column_rank([[positions[k] for k in col] for col in columns])
    e = len(identities)
    return {"dim": e + tau - 2 * r, "e": e, "tau": tau, "r": r, "factors": factors}


FACTORS = set()


def pair_dim(earring, blue):
    result = mapping_cone(earring, blue)
    FACTORS.update(result["factors"])
    return result["dim"]


def smoothing(translated, lifts, pairing):
    p2, mismatch = iota_pairing(edges_of(translated), lifts[0], pairing, lifts[1], EPS)
    if mismatch > 1.0e-8:
        raise AssertionError("smoothing is not iota-invariant")
    components = quiet(smooth_curve, translated,
                       [(lifts[0], pairing), (lifts[1], p2)], EPS)
    return components


def main():
    start = time.time()
    _, blue, xinfo = quiet(q7.build_q7)
    rows = census.provenance_census(blue, xinfo)
    report("the model has 82 self-intersection orbits", len(rows) == 82)
    translated = q7.translate_to_kwz_chart(blue)
    data = q7.encode_type_d(translated)
    data["precurve"] = q7.geometric_precurve(data)
    base = (data, data["precurve"], data["precurve"]["delta"])
    by_signature = {census.orbit_signature(pre): (pt, pre)
                    for pt, pre in orbit_group(translated)}

    earrings = {label: quiet(closure.rational_earring_type_d, slope)
                for label, slope in (("E", Fraction(-1, 2)),
                                     ("E*", Fraction(-3, 4)),
                                     ("E3/4", Fraction(3, 4)))}
    print("earring generators: " + ", ".join(
        f"{label} {len(pre['vertices'])}" for label, (_, pre) in earrings.items()))
    n_dims = tuple(pair_dim(earrings[x], base) for x in ("E", "E*", "E3/4"))
    print(f"N: 31 generators; dimensions against (E, E*, E3/4) = {n_dims}")
    report("N pairs to (7,25) with E and E*", n_dims[:2] == (7, 25))
    reference = tuple(
        q7.exact_red_blue_homology(red[0], red[1]["vertices"], red[1]["delta"],
                                   data, data["precurve"], data["precurve"]["delta"])["homology"]
        for red in (earrings["E"], earrings["E*"], earrings["E3/4"]))
    report("the reduction agrees with q7_kwz.exact_red_blue_homology on N", reference == n_dims)

    census_pairs, opposite_pairs, excluded = {}, {}, []
    opposite_data = {}
    for row in rows:
        support = by_signature[row["signature"]]
        lifts = support[1]
        switch = quiet(census.physical_switch_only, translated, support, data, 1)
        if not switch["module_matches"]:
            excluded.append(row["name"])
        for pairing, store in ((1, census_pairs), (2, opposite_pairs)):
            obj = standalone(smoothing(translated, lifts, pairing))
            dims = (pair_dim(earrings["E"], obj), pair_dim(earrings["E*"], obj))
            store[row["name"]] = (len(obj[1]["vertices"]), dims)
            if pairing == 2:
                opposite_data[row["name"]] = obj

    print("\ncensus smoothings that change the module of N:")
    for name in excluded:
        gens, dims = census_pairs[name]
        print(f"  {name}: {gens} generators; dimensions (E, E*) = {dims}")
    report("exactly the nine orbits S7,S8,S23,S24,S49,S52,S53,S61,S73 change the module",
           set(excluded) == set(EXCLUDED) and len(excluded) == 9)
    report("each of the nine encodings has 31 generators",
           all(census_pairs[n][0] == 31 for n in EXCLUDED))
    report("the nine dimension pairs are (3,23) at S7, (5,23) at S8,S53,S73, (7,25) otherwise",
           all(census_pairs[n][1] == EXCLUDED[n] for n in EXCLUDED))

    census_hist = Counter(dims for _, dims in census_pairs.values())
    print(f"\nall 82 census smoothings, (E, E*) histogram: {dict(sorted(census_hist.items()))}")
    census_hits = [n for n, (_, d) in census_pairs.items() if d == (9, 31)]
    report("among the 82 census smoothings, (9,31) occurs only at S18 and S25",
           census_hits == SELECTED)
    named = {n: census_pairs[n][1] for n in ("S18", "S25", "S69", "S74")}
    report("the named census smoothings give (9,31),(9,31),(9,23),(9,25)",
           named == {"S18": (9, 31), "S25": (9, 31), "S69": (9, 23), "S74": (9, 25)})

    opposite_hist = Counter(dims for _, dims in opposite_pairs.values())
    print(f"\nall 82 opposite smoothings, (E, E*) histogram: {dict(sorted(opposite_hist.items()))}")
    opposite_hits = [n for n, (_, d) in opposite_pairs.items() if d == (9, 31)]
    report("among the 82 opposite smoothings, (9,31) occurs only at S18 and S25",
           opposite_hits == SELECTED)

    print("\nopposite smoothing at S18 and S25 against the slope-3/4 earring:")
    s18 = q7.collect_f2(data["precurve"]["delta"] + quiet(
        census.physical_switch_only, translated,
        by_signature[next(r["signature"] for r in rows if r["name"] == "S18")],
        data, 1)["switch"])
    n18 = pair_dim(earrings["E3/4"], (data, data["precurve"], s18))
    print(f"  (N, b_S18): {n18}")
    for name in SELECTED:
        obj = opposite_data[name]
        dim34 = pair_dim(earrings["E3/4"], obj)
        print(f"  {name}: pillowcase components {obj[3]}, "
              f"{len(obj[1]['vertices'])} generators, E3/4 dimension {dim34}")
        report(f"opposite smoothing at {name}: one component, 21 generators, dimension 69",
               obj[3] == 1 and len(obj[1]["vertices"]) == 21 and dim34 == 69)
    report("(N, b_S18) has dimension 103 against the slope-3/4 earring", n18 == 103)

    print(f"\ninvariant factors met in all pairings: {sorted(bin(f)[2:] for f in FACTORS)}"
          " (binary, lowest power of U on the right)")
    report("every invariant factor of every face complex equals U", FACTORS == {2})

    slopes = (Fraction(3, 4), Fraction(1, 3), Fraction(1, 7))
    determinant = abs(4 * 3 * 7 * sum(slopes))
    report("num(Q_{3/4}+Q_{1/3}+Q_{1/7}) has determinant 103", determinant == 103)

    print(f"\nelapsed {time.time() - start:.1f} s")
    if FAILURES:
        print(f"{len(FAILURES)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
