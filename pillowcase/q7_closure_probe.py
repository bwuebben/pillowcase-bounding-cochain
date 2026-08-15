#!/usr/bin/env python3
"""Exact rational-earring probes for the three q=7 deformation classes.

The fixed blue object models the perturbed Conway sum Q_{1/3}+Q_{1/7}.
Pair it with the earring of a rational tangle Q_r.  For r=1/n the corresponding
numerator closure is the pretzel knot P(n,3,7).  The slope r=-3/4 gives the
auxiliary Montesinos knot M(-3/4,1/3,1/7), which distinguishes the three
displayed deformation classes.

This script is an algebraic discriminator only.  Interpreting a row as
instanton homology still requires the CHKK tangle-pairing conjecture.
"""

from __future__ import annotations

import argparse
import contextlib
import io
from fractions import Fraction

from bigons import simplify
from earring import f8
import q7_kwz as q7


AUX_SLOPE = Fraction(-3, 4)
AUX_PD = [
    (24, 13, 25, 14), (10, 25, 11, 26), (26, 11, 27, 12),
    (12, 27, 13, 0), (7, 14, 8, 15), (15, 8, 16, 9),
    (9, 16, 10, 17), (23, 6, 24, 7), (5, 22, 6, 23),
    (21, 4, 22, 5), (3, 20, 4, 21), (19, 2, 20, 3),
    (1, 18, 2, 19), (17, 0, 18, 1),
]
AUX_KHOCA_GAUSS = [18, 16, 11, 24, 10, 23, 9, 22, 14, 21, 20, 19]
# (homological degree, quantum degree, multiplicity), as printed by Khoca 1.5
# over F_2.  Both AUX_PD and AUX_KHOCA_GAUSS give the same total rank.
AUX_REDUCED_KH = [
    (0, -12, 1), (2, -16, 1), (3, -18, 1), (4, -18, 1),
    (5, -22, 1), (5, -20, 1), (6, -22, 3), (7, -24, 3),
    (8, -26, 4), (9, -28, 5), (10, -30, 4), (11, -32, 3),
    (12, -34, 2), (13, -36, 1),
]
# Spherogram 2.4.1 Seifert matrix for AUX_PD.
AUX_SEIFERT_MATRIX = [
    [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1],
]
# Coefficients of det(V-tV^T), in ascending powers of t.
AUX_ALEXANDER = (1, -1, 0, 2, -4, 5, -5, 5, -4, 2, 0, -1, 1)
EXPECTED_PHYSICAL_SPAN_RANKS = {
    (): (7, 25),
    ("S18",): (9, 31),
    ("S25",): (9, 31),
    ("S18", "S25"): (7, 25),
    ("S69",): (9, 23),
    ("S18", "S69"): (7, 25),
    ("S25", "S69"): (7, 25),
    ("S18", "S25", "S69"): (9, 23),
    ("S74",): (9, 25),
    ("S18", "S74"): (9, 25),
    ("S25", "S74"): (7, 25),
    ("S18", "S25", "S74"): (7, 25),
    ("S69", "S74"): (7, 23),
    ("S18", "S69", "S74"): (7, 23),
    ("S25", "S69", "S74"): (9, 23),
    ("S18", "S25", "S69", "S74"): (9, 23),
}


def integer_determinant(matrix):
    """Exact Bareiss determinant of an integer matrix."""
    matrix = [row[:] for row in matrix]
    size = len(matrix)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if not matrix[pivot_index][pivot_index]:
            swap = next(
                row for row in range(pivot_index + 1, size)
                if matrix[row][pivot_index]
            )
            matrix[pivot_index], matrix[swap] = (
                matrix[swap], matrix[pivot_index])
            sign = -sign
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index]
                    * matrix[pivot_index][column]
                )
                matrix[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def auxiliary_alexander_certificate():
    """Verify det(V-tV^T) and return its norm and determinant."""
    size = len(AUX_SEIFERT_MATRIX)
    for value in range(-6, 7):
        matrix = [[
            AUX_SEIFERT_MATRIX[row][column]
            - value * AUX_SEIFERT_MATRIX[column][row]
            for column in range(size)
        ] for row in range(size)]
        determinant = integer_determinant(matrix)
        polynomial = sum(
            coefficient * value ** degree
            for degree, coefficient in enumerate(AUX_ALEXANDER)
        )
        if determinant != polynomial:
            raise AssertionError("auxiliary Alexander determinant mismatch")
    return {
        "coefficients": AUX_ALEXANDER,
        "norm": sum(abs(value) for value in AUX_ALEXANDER),
        "determinant": abs(sum(
            coefficient * (-1) ** degree
            for degree, coefficient in enumerate(AUX_ALEXANDER)
        )),
    }


def verify_external_packages():
    """Optionally rerun the PD, Seifert-matrix, and Khoca certificates."""
    from khoca import InteractiveCalculator
    from spherogram import RationalTangle

    tangle = (
        RationalTangle(-3, 4)
        + RationalTangle(1, 3)
        + RationalTangle(1, 7)
    )
    link = tangle.numerator_closure()
    if link.PD_code() != AUX_PD:
        raise AssertionError("Spherogram auxiliary PD changed")
    if link.seifert_matrix() != AUX_SEIFERT_MATRIX:
        raise AssertionError("Spherogram auxiliary Seifert matrix changed")
    result, _ = InteractiveCalculator(2)(AUX_PD, print_messages=True)
    if sum(term[3] for term in result[0]) != 31:
        raise AssertionError("Khoca auxiliary reduced rank is not 31")


def rational_earring_type_d(slope: Fraction):
    """Return the canceled type-D object for the mirrored Q_slope earring."""
    if not slope:
        direction = (1, 0)
    else:
        direction = (slope.denominator, -slope.numerator)
    red = simplify(
        f8(direction, eps=q7.DEFAULT["red_eps"],
           phi=q7.DEFAULT["red_phi"])[0],
        2e-4,
    )
    red = q7.translate_to_kwz_chart(red)
    data = q7.encode_type_d(red)
    pre = q7.cyclic_precurve(data)
    if pre["residue"]:
        raise AssertionError(
            f"Q_{{{slope}}} red differential does not square to zero")
    return data, pre


def q7_deformation_objects():
    """Build the common module and one differential for each displayed class."""
    with contextlib.redirect_stdout(io.StringIO()):
        red, blue, _ = q7.build_q7()
        red = q7.translate_to_kwz_chart(red)
        blue = q7.translate_to_kwz_chart(blue)
        targets = {
            name: q7.translate_target_to_kwz(point)
            for name, point in q7.TARGETS.items()
        }
        data = q7.report_encoding(blue, targets["S69"])
        certificate = q7.report_admissible_kwz(red, blue, targets, data)
    return data, {
        "undeformed": data["precurve"]["delta"],
        "S18/S25": certificate["physical"]["S18"]["delta"],
        "S69": certificate["physical"]["S69"]["delta"],
        "S74": certificate["physical"]["S74"]["delta"],
    }


def q7_physical_switches():
    """Return the four exact physical smoothing switches on the common module."""
    with contextlib.redirect_stdout(io.StringIO()):
        red, blue, _ = q7.build_q7()
        red = q7.translate_to_kwz_chart(red)
        blue = q7.translate_to_kwz_chart(blue)
        targets = {
            name: q7.translate_target_to_kwz(point)
            for name, point in q7.TARGETS.items()
        }
        data = q7.report_encoding(blue, targets["S69"])
        certificate = q7.report_admissible_kwz(red, blue, targets, data)
    return data, {
        name: result["b"]
        for name, result in certificate["physical"].items()
    }


def physical_switch_span_census():
    """Test every F2-sum of the four physical smoothing switches exactly."""
    blue_data, switches = q7_physical_switches()
    blue_pre = blue_data["precurve"]
    names = tuple(sorted(switches))
    earrings = {}
    for label, slope in (("original", Fraction(-1, 2)),
                         ("auxiliary", AUX_SLOPE)):
        earrings[label] = rational_earring_type_d(slope)

    rows = []
    for mask in range(1 << len(names)):
        selected = tuple(
            name for index, name in enumerate(names) if mask & (1 << index)
        )
        switch = q7.collect_f2([
            term for name in selected for term in switches[name]
        ])
        residue = q7.mc_residue(blue_pre["delta"], switch)[0]
        row = {
            "mask": mask,
            "selected": selected,
            "terms": len(switch),
            "maurer_cartan": not residue,
            "residue": residue,
            "ranks": {},
            "certificates": {},
        }
        if not residue:
            blue_delta = q7.collect_f2(blue_pre["delta"] + switch)
            for label, (red_data, red_pre) in earrings.items():
                result = q7.exact_red_blue_homology(
                    red_data, red_pre["vertices"], red_pre["delta"],
                    blue_data, blue_pre, blue_delta,
                )
                row["ranks"][label] = result["homology"]
                row["certificates"][label] = tuple(
                    result[key] for key in (
                        "identities", "tail_homology",
                        "beta_rank", "homology",
                    )
                )
        rows.append(row)
    return rows


def earring_mapping_cone():
    """Construct CHKK's [I->I](N)=Cone(H id_N) in three-face notation.

    On a left-idempotent generator the central element H is L+M^2; on a
    right-idempotent generator it is R+M^2.  The returned 62-generator type-D
    object is purely algebraic and does not assume CHKK Conjecture F.
    """
    blue_data, _ = q7_deformation_objects()
    blue_pre = blue_data["precurve"]
    vertices = tuple(blue_pre["vertices"])
    offset = max(vertices) + 1
    second = tuple(vertex + offset for vertex in vertices)
    unique = [dict(record) for record in blue_data["unique"]]
    unique.extend(dict(record) for record in blue_data["unique"])
    delta = []
    for copy_offset in (0, offset):
        delta.extend(
            (source + copy_offset, word, target + copy_offset)
            for source, word, target in blue_pre["delta"]
        )
    for vertex in vertices:
        arc = blue_data["unique"][vertex]["arc"]
        outside_face = "L" if arc == "L" else "R"
        delta.extend([
            (vertex, (outside_face, 1), vertex + offset),
            (vertex, ("M", 2), vertex + offset),
        ])
    delta = q7.collect_f2(delta)
    residue = q7.apply_mul(delta, delta)
    if residue:
        raise AssertionError("CHKK earring mapping cone does not square to zero")
    data = dict(blue_data)
    data["unique"] = unique
    pre = {
        "vertices": vertices + second,
        "delta": delta,
        "residue": residue,
    }
    return data, pre


def corrected_end_cohomology_profile():
    """Compute all quantum sectors of corrected-chart H^(q,-1) End(N)."""
    data, _ = q7_deformation_objects()
    pre = data["precurve"]
    grades, consistent = q7.generator_bigradings(
        pre["vertices"], pre["delta"])
    if not consistent:
        raise AssertionError("corrected q=7 generator grading is inconsistent")
    degree_minus_one = q7.end_terms_in_delta_degree(data, pre, grades, -1)
    quantum_degrees = sorted({
        q7.morphism_term_bigrading(term, grades)[0]
        for term in degree_minus_one
    })
    sectors = {}
    for quantum_degree in quantum_degrees:
        result = q7.end_cohomology_data(data, quantum_degree)
        if result["representatives"]:
            sectors[quantum_degree] = result
    return sectors


def corrected_representative_mc_profile():
    """Square-zero census in the dynamically computed 20-class H^-1 basis."""
    sectors = corrected_end_cohomology_profile()
    labels = []
    representatives = []
    for quantum_degree, result in sectors.items():
        for index, representative in enumerate(result["representatives"]):
            labels.append((quantum_degree, index))
            representatives.append(representative)
    masks = q7._square_zero_masks(representatives)
    return {
        "labels": tuple(labels),
        "representatives": tuple(representatives),
        "masks": tuple(masks),
    }


def physical_switch_cohomology_coordinates():
    """Coordinates of the four switches in the corrected H^-1 End basis."""
    data, switches = q7_physical_switches()
    sectors = corrected_end_cohomology_profile()
    grades = next(iter(sectors.values()))["grades"]
    coordinates = {}
    for name, switch in switches.items():
        by_quantum = {}
        for term in switch:
            degree = q7.morphism_term_bigrading(term, grades)[0]
            by_quantum.setdefault(degree, []).append(term)
        coordinate = []
        for degree, component in sorted(by_quantum.items()):
            result = sectors.get(degree)
            if result is None:
                result = q7.end_cohomology_data(data, degree)
            basis = result["bases"][-1]
            index = {term: position for position, term in enumerate(basis)}
            component_bits = sum(1 << index[term] for term in component)
            boundary_values = []
            for term in result["bases"][0]:
                image = q7.end_differential(data["precurve"]["delta"], [term])
                boundary_values.append(sum(1 << index[value] for value in image))
            boundaries = q7._bit_span_basis(boundary_values)
            representative_bits = [
                sum(1 << index[term] for term in representative)
                for representative in result["representatives"]
            ]
            matches = []
            for mask in range(1 << len(representative_bits)):
                value = component_bits
                for position, representative in enumerate(representative_bits):
                    if mask & (1 << position):
                        value ^= representative
                if q7._bit_in_span(value, boundaries):
                    matches.append(mask)
            if len(matches) != 1:
                raise AssertionError(f"{name} cohomology coordinates are ambiguous")
            mask = matches[0]
            coordinate.extend(
                (degree, position)
                for position in range(len(representative_bits))
                if mask & (1 << position)
            )
        coordinates[name] = tuple(coordinate)
    return coordinates


def probe(slopes):
    blue_data, objects = q7_deformation_objects()
    blue_pre = blue_data["precurve"]
    rows = []
    for slope in slopes:
        red_data, red_pre = rational_earring_type_d(slope)
        ranks = {}
        certificates = {}
        for name, delta in objects.items():
            result = q7.exact_red_blue_homology(
                red_data, red_pre["vertices"], red_pre["delta"],
                blue_data, blue_pre, delta,
            )
            ranks[name] = result["homology"]
            certificates[name] = tuple(result[key] for key in (
                "identities", "tail_homology", "beta_rank", "homology"
            ))
        rows.append({
            "slope": slope,
            "determinant": abs(
                21 * slope.numerator + 10 * slope.denominator),
            "red_vertices": len(red_pre["vertices"]),
            "ranks": ranks,
            "certificates": certificates,
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", nargs="*", type=int,
                        help="nonzero pretzel parameters (default: -5..8)")
    parser.add_argument(
        "--slope", action="append", type=Fraction, default=[],
        help="general rational tangle slope p/q; may be repeated",
    )
    parser.add_argument(
        "--selection-certificate", action="store_true",
        help="verify the r=-3/4 auxiliary-closure selection certificate",
    )
    parser.add_argument(
        "--external", action="store_true",
        help="also rerun Spherogram 2.4.1 and Khoca 1.5",
    )
    parser.add_argument(
        "--physical-span-census", action="store_true",
        help="test all sixteen sums of the four physical q=7 switches",
    )
    parser.add_argument(
        "--earring-cone", action="store_true",
        help="construct CHKK's 62-generator [I->I](N) mapping cone",
    )
    parser.add_argument(
        "--corrected-end-profile", action="store_true",
        help="recompute every corrected-chart H^(q,-1) End(N) sector",
    )
    parser.add_argument(
        "--corrected-mc-profile", action="store_true",
        help="test square-zero sums in the corrected 20-class End basis",
    )
    parser.add_argument(
        "--physical-cohomology", action="store_true",
        help="locate the four physical switches in corrected H^-1 End(N)",
    )
    args = parser.parse_args()
    if any(n == 0 for n in args.n):
        parser.error("n must be nonzero")
    slopes = args.slope or [Fraction(1, n) for n in (
        args.n or [n for n in range(-5, 9) if n]
    )]

    rows = probe(slopes)
    print("slope   det   |red|   undeformed   S18/S25   S69   S74")
    for row in rows:
        ranks = row["ranks"]
        print(f"{str(row['slope']):>5}  {row['determinant']:>4}"
              f"  {row['red_vertices']:>5}"
              f"  {ranks['undeformed']:>11}  {ranks['S18/S25']:>8}"
              f"  {ranks['S69']:>4}  {ranks['S74']:>4}")

    if args.selection_certificate:
        row = next((row for row in rows if row["slope"] == AUX_SLOPE), None)
        if row is None:
            row = probe([AUX_SLOPE])[0]
        expected = {
            "undeformed": 25, "S18/S25": 31, "S69": 23, "S74": 25,
        }
        if row["ranks"] != expected:
            raise AssertionError("auxiliary closure pairing ranks changed")
        alexander = auxiliary_alexander_certificate()
        kh_rank = sum(multiplicity for _, _, multiplicity in AUX_REDUCED_KH)
        if alexander["norm"] != kh_rank or kh_rank != 31:
            raise AssertionError("auxiliary instanton squeeze changed")
        span = physical_switch_span_census()
        span_ranks = {
            candidate["selected"]: (
                candidate["ranks"].get("original"),
                candidate["ranks"].get("auxiliary"),
            )
            for candidate in span
        }
        if (not all(candidate["maurer_cartan"] for candidate in span)
                or span_ranks != EXPECTED_PHYSICAL_SPAN_RANKS):
            raise AssertionError("physical-switch span certificate changed")
        if args.external:
            verify_external_packages()
        print("selection certificate: PASS")
        print(f"  exact ranks={row['ranks']}")
        print(f"  identity/tail/beta/homology={row['certificates']}")
        print(f"  Alexander norm={alexander['norm']}; "
              f"determinant={alexander['determinant']}")
        print(f"  reduced Kh(F2) rank={kh_rank}")
        print("  physical span=16 Maurer-Cartan sums; "
              "rank pair (9,31) only for S18 and S25")
        print("  conditional survivor=S18/S25")

    if args.physical_span_census:
        print("physical-switch span over F2:")
        print("  subset                         |b|   MC   original   auxiliary")
        for row in physical_switch_span_census():
            subset = "+".join(row["selected"]) or "0"
            ranks = row["ranks"]
            original = str(ranks.get("original", "-"))
            auxiliary = str(ranks.get("auxiliary", "-"))
            print(f"  {subset:<30} {row['terms']:>3}  "
                  f"{'yes' if row['maurer_cartan'] else 'no ':>3}"
                  f"  {original:>9}  {auxiliary:>10}")

    if args.earring_cone:
        _, cone = earring_mapping_cone()
        print("CHKK earring mapping cone:")
        print(f"  vertices={len(cone['vertices'])}; "
              f"arrows={len(cone['delta'])}; residue={len(cone['residue'])}")

    if args.corrected_end_profile:
        sectors = corrected_end_cohomology_profile()
        dimensions = {
            int(degree): len(result["representatives"])
            for degree, result in sectors.items()
        }
        print("corrected-chart degree-minus-one End cohomology:")
        print(f"  dimensions={dimensions}")
        print(f"  total dimension={sum(dimensions.values())}")

    if args.corrected_mc_profile:
        profile = corrected_representative_mc_profile()
        total = 1 << len(profile["representatives"])
        print("corrected-chart cohomology-representative MC census:")
        print(f"  square-zero sums={len(profile['masks'])} of {total}")
        print(f"  nonzero square sums={total - len(profile['masks'])}")

    if args.physical_cohomology:
        print("physical-switch cohomology coordinates:")
        for name, coordinate in physical_switch_cohomology_coordinates().items():
            print(f"  {name}: {coordinate}")


if __name__ == "__main__":
    main()
