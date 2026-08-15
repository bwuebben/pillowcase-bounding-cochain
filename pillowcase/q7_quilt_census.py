#!/usr/bin/env python3
"""Finite Gate-1 census for the q=7 Conway-sum object.

This script records the part of a proposed quilt/Yoneda calculation that is
already finite geometry.  It does *not* count holomorphic quilts and it does
not prove the CHKK instanton--pillowcase correspondence.

The certificate has two layers.

1. Every self-intersection of the resolved curve for Q_{1/3}+Q_{1/7} is
   classified by provenance: main--main, connector--main, or
   connector--connector.  Connector branches retain the seam-fibre circle,
   U/D lift, and connector number from ``resolve.py``.

2. Pairing-1 surgery at every node orbit is encoded in the admissible
   (0,pi)-special KWZ chart.  We record whether it preserves the cancelled
   31-generator module and, if so, whether its differential change is a
   Maurer--Cartan element.  The connector--main subcensus is retained for the
   unrestricted-span tests.  This is an exhaustive census of physical
   reconnections, not a theorem that figure-eight bubbles have singleton
   support.

The missing analytic implication is deliberately printed at the end:

    rigid C_3 figure-eight output
        => one module-preserving physical reconnection (or the named span).

Run from ``pillowcase`` or from the repository root:

    python3 pillowcase/q7_quilt_census.py
    python3 pillowcase/q7_quilt_census.py --switch-census
    python3 pillowcase/q7_quilt_census.py --pairing-census
    python3 pillowcase/q7_quilt_census.py --all-node-census
    python3 pillowcase/q7_quilt_census.py --two-switch-census
    python3 pillowcase/q7_quilt_census.py --two-switch-census \
        --strict-pair-census
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import q7_kwz as q7
import q7_closure_probe as closure
from bigons import _tdist, edges_of
from earring import P_point
from maurer_cartan import orbit_group
from resolve import self_crossings_T2
from surgery_check import iota_pairing, smooth_curve
from tangles import PI, TAU


DISTANCE_TOL = 2.0e-3
ANGLE_TOL = math.radians(2.0)

EXPECTED_PROVENANCE = {
    "connector--connector": 30,
    "connector--main": 52,
}
EXPECTED_CANONICAL = {
    0: "S19",
    1: "S31",
    2: "S45",
    3: "S80",
    4: "S8",
    5: "S53",
}
EXPECTED_MODULE_EXCLUSIONS = {"S7", "S23", "S24", "S49", "S52", "S61"}
EXPECTED_ALL_MODULE_EXCLUSIONS = EXPECTED_MODULE_EXCLUSIONS | {"S8", "S53", "S73"}
EXPECTED_DUPLICATE_SWITCHES = {
    frozenset(pair) for pair in (
        ("S0", "S14"),
        ("S5", "S11"),
        ("S12", "S16"),
        ("S26", "S33"),
        ("S46", "S48"),
        ("S50", "S58"),
        ("S54", "S57"),
        ("S70", "S72"),
    )
}
EXPECTED_NODE_RANK_CENSUS = Counter({
    (5, 23): 1,
    (5, 25): 1,
    (7, 23): 5,
    (7, 25): 25,
    (7, 27): 5,
    (9, 23): 1,
    (9, 25): 6,
    (9, 31): 2,
})
EXPECTED_DISTINCT_RANK_CENSUS = Counter({
    (5, 23): 1,
    (5, 25): 1,
    (7, 23): 4,
    (7, 25): 19,
    (7, 27): 4,
    (9, 23): 1,
    (9, 25): 6,
    (9, 31): 2,
})
EXPECTED_ALL_NODE_RANK_CENSUS = Counter({
    (5, 23): 1,
    (5, 25): 2,
    (7, 23): 12,
    (7, 25): 39,
    (7, 27): 8,
    (9, 23): 1,
    (9, 25): 8,
    (9, 31): 2,
})
EXPECTED_ALL_DISTINCT_RANK_CENSUS = Counter({
    (5, 23): 1,
    (5, 25): 2,
    (7, 23): 6,
    (7, 25): 23,
    (7, 27): 4,
    (9, 23): 1,
    (9, 25): 6,
    (9, 31): 2,
})
EXPECTED_ALL_DUPLICATE_SWITCHES = {
    frozenset(names) for names in (
        ("S0", "S14", "S15", "S21"),
        ("S1", "S2"),
        ("S5", "S6", "S11"),
        ("S12", "S16", "S17"),
        ("S26", "S27", "S28", "S32", "S33"),
        ("S36", "S37"),
        ("S41", "S42", "S44"),
        ("S46", "S47", "S48"),
        ("S50", "S51", "S56", "S58"),
        ("S54", "S55", "S57"),
        ("S63", "S64", "S65"),
        ("S70", "S71", "S72"),
        ("S77", "S78", "S79"),
    )
}
EXPECTED_TWO_SWITCH_RANK_CENSUS = Counter({
    (5, 23): 35,
    (5, 25): 18,
    (7, 23): 118,
    (7, 25): 205,
    (7, 27): 78,
    (9, 23): 32,
    (9, 25): 126,
    (9, 27): 35,
    (9, 29): 5,
    (9, 31): 41,
    (11, 25): 2,
    (11, 27): 4,
    (11, 31): 2,
    (11, 33): 2,
})


def _short_delta(a, b):
    """Shortest coordinate difference from ``a`` to ``b`` on T^2."""
    return ((b[0] - a[0] + PI) % TAU - PI,
            (b[1] - a[1] + PI) % TAU - PI)


def _point_segment_distance(point, start, end):
    """Distance from a torus point to one short torus segment."""
    direction = _short_delta(start, end)
    offset = _short_delta(start, point)
    norm2 = direction[0] ** 2 + direction[1] ** 2
    if norm2 == 0.0:
        return math.hypot(*offset)
    parameter = max(0.0, min(1.0,
        (offset[0] * direction[0] + offset[1] * direction[1]) / norm2))
    residual = (offset[0] - parameter * direction[0],
                offset[1] - parameter * direction[1])
    return math.hypot(*residual)


def _unoriented_angle(first, second):
    """Angle between unoriented lines, in [0,pi/2]."""
    n1 = math.hypot(*first)
    n2 = math.hypot(*second)
    if n1 == 0.0 or n2 == 0.0:
        return PI / 2
    cosine = abs((first[0] * second[0] + first[1] * second[1]) / (n1 * n2))
    return math.acos(max(-1.0, min(1.0, cosine)))


def connector_records(xinfo):
    """Flatten ``resolve`` connector records while retaining their provenance."""
    records = []
    for entry in xinfo:
        for lift in ("U", "D"):
            for number, polyline in enumerate(entry[lift]):
                records.append({
                    "circle": entry["circle"],
                    "seam": entry["seam"],
                    "lift": lift,
                    "number": number,
                    "polyline": polyline,
                })
    return records


def classify_branch(point, tangent, connectors):
    """Return the connector aligned with a crossing branch, or ``None``.

    Spatial distance alone cannot distinguish the two branches at a
    connector--strand crossing: both pass through the same point.  The
    connector is therefore selected by its tangent line as well.
    """
    candidates = []
    for connector in connectors:
        polyline = connector["polyline"]
        for start, end in zip(polyline, polyline[1:]):
            distance = _point_segment_distance(point, start, end)
            if distance > DISTANCE_TOL:
                continue
            direction = _short_delta(start, end)
            angle = _unoriented_angle(tangent, direction)
            candidates.append((angle, distance, connector))
    if not candidates:
        return None
    angle, distance, connector = min(candidates, key=lambda item: item[:2])
    if angle > ANGLE_TOL:
        return None
    return {
        "circle": connector["circle"],
        "seam": connector["seam"],
        "lift": connector["lift"],
        "number": connector["number"],
        "distance": distance,
        "angle": angle,
    }


def orbit_signature(preimages):
    """Coordinate-independent edge signature of a pillowcase crossing orbit."""
    return tuple(sorted(
        (min(crossing["kA"], crossing["kB"]),
         max(crossing["kA"], crossing["kB"]))
        for crossing in preimages
    ))


def provenance_census(blue, xinfo):
    connectors = connector_records(xinfo)
    rows = []
    for index, (point, preimages) in enumerate(orbit_group(blue)):
        lifts = []
        for crossing in preimages:
            branch_a = classify_branch(crossing["pt"], crossing["dirA"], connectors)
            branch_b = classify_branch(crossing["pt"], crossing["dirB"], connectors)
            count = int(branch_a is not None) + int(branch_b is not None)
            kind = ("main--main", "connector--main", "connector--connector")[count]
            lifts.append({
                "kind": kind,
                "A": branch_a,
                "B": branch_b,
                "edges": tuple(sorted((crossing["kA"], crossing["kB"]))),
            })
        kinds = {lift["kind"] for lift in lifts}
        rows.append({
            "name": f"S{index}",
            "index": index,
            "point": tuple(point),
            "signature": orbit_signature(preimages),
            "kind": next(iter(kinds)) if len(kinds) == 1 else "inconsistent",
            "lifts": lifts,
            "preimages": preimages,
        })
    return rows


def locate_named(rows):
    named = {}
    for name, target in q7.TARGETS.items():
        row = min(rows, key=lambda item: _tdist(item["point"], target))
        if _tdist(row["point"], target) > 0.02:
            raise AssertionError(f"could not locate {name}")
        named[name] = row
    return named


def canonical_connector_orbits(rows, xinfo):
    """Match the intrinsic connector--connector node of each fibre circle."""
    result = {}
    for entry in xinfo:
        points = self_crossings_T2([entry["U"][0], entry["U"][1]])
        if len(points) != 1:
            raise AssertionError(
                f"circle {entry['circle']} has {len(points)} U connector crossings")
        pillow_point = P_point(points[0])
        row = min(rows, key=lambda item: _tdist(item["point"], pillow_point))
        if _tdist(row["point"], pillow_point) > 3.0e-3:
            raise AssertionError(
                f"canonical crossing of circle {entry['circle']} was not found")
        result[entry["circle"]] = row["name"]
    return result


def _component_signature(component_data):
    return tuple(sorted(
        (dot["arc"], round(dot["y"], 7))
        for dot in component_data["unique"]
    ))


def physical_switch_only(blue, support, data, pairing=1, eps=0.006):
    """Encode a physical smoothing without doing any closure-Hom calculation."""
    point, lifts = support
    if len(lifts) != 2:
        return {"module_matches": False, "error": f"{len(lifts)} torus lifts"}
    p2, mismatch = iota_pairing(edges_of(blue), lifts[0], pairing,
                                lifts[1], eps)
    if mismatch > 1.0e-8:
        return {"module_matches": False, "error": f"iota mismatch {mismatch}"}
    components = smooth_curve(blue, [(lifts[0], pairing), (lifts[1], p2)], eps)

    encoded = []
    quotient_groups = defaultdict(list)
    for index, component in enumerate(components):
        component_data = q7.encode_type_d(component)
        invariant = len(component_data["hits"]) != len(component_data["unique"])
        component_pre = (q7.geometric_precurve(component_data)
                         if invariant else q7.cyclic_precurve(component_data))
        signature = _component_signature(component_data)
        encoded.append((component_data, component_pre, invariant, signature))
        quotient_groups[(invariant, signature)].append(index)

    selected = []
    for (invariant, _), indices in quotient_groups.items():
        if invariant:
            selected.extend(indices)
        elif len(indices) == 2:
            selected.append(indices[0])
        else:
            return {
                "module_matches": False,
                "error": "non-invariant components do not form an iota pair",
            }

    mapped_vertices = []
    mapped_delta = []
    for index in selected:
        component_data, component_pre, _, _ = encoded[index]
        try:
            matching = q7._match_component_dots(component_data, data)
        except AssertionError as error:
            return {"module_matches": False, "error": str(error)}
        mapped_vertices.extend(matching[i] for i in component_pre["vertices"])
        mapped_delta.extend(
            (matching[i], word, matching[j])
            for i, word, j in component_pre["delta"]
        )

    base = data["precurve"]
    module_matches = (
        len(mapped_vertices) == len(set(mapped_vertices))
        and set(mapped_vertices) == set(base["vertices"])
    )
    if not module_matches:
        return {"module_matches": False, "error": "31-dot module changed"}
    mapped_delta = q7.collect_f2(mapped_delta)
    switch = q7.collect_f2(base["delta"] + mapped_delta)
    residue, differential, square = q7.mc_residue(base["delta"], switch)
    return {
        "module_matches": True,
        "pairings": (pairing, p2),
        "terms": len(switch),
        "switch": switch,
        "mc": not residue,
        "d_terms": len(differential),
        "square_terms": len(square),
        "residue_terms": len(residue),
    }


def switch_census(rows, blue, pairing_ranks=False, all_nodes=False):
    """Encode singleton smoothings.

    By default this preserves the original connector--main subcensus used by
    the two-switch span calculation.  ``all_nodes=True`` tests all 82 node
    orbits, including connector--connector nodes.
    """
    translated = q7.translate_to_kwz_chart(blue)
    data = q7.encode_type_d(translated)
    data["precurve"] = q7.geometric_precurve(data)
    if (len(data["precurve"]["vertices"]), len(data["precurve"]["delta"])) != (31, 30):
        raise AssertionError("admissible q=7 base object is not the expected 31/30 complex")
    grades, grading_ok = q7.generator_bigradings(
        data["precurve"]["vertices"], data["precurve"]["delta"])
    if not grading_ok:
        raise AssertionError("admissible q=7 module has no relative KWZ grading")

    translated_by_signature = {
        orbit_signature(preimages): (point, preimages)
        for point, preimages in orbit_group(translated)
    }
    earrings = {}
    if pairing_ranks:
        for label, slope in (("original", Fraction(-1, 2)),
                             ("auxiliary", closure.AUX_SLOPE)):
            earrings[label] = closure.rational_earring_type_d(slope)
    results = []
    for row in rows:
        if not all_nodes and row["kind"] != "connector--main":
            continue
        support = translated_by_signature.get(row["signature"])
        if support is None:
            results.append((row, {"module_matches": False,
                                  "error": "translated orbit not found"}))
            continue
        try:
            result = physical_switch_only(translated, support, data)
        except (AssertionError, KeyError, ValueError) as error:
            result = {"module_matches": False, "error": str(error)}
        if result["module_matches"]:
            result["degrees"] = [
                q7.morphism_term_bigrading(term, grades)
                for term in result["switch"]
            ]
            result["delta_homogeneous"] = {
                degree[1] for degree in result["degrees"]
            } == {-1}
            result["short_decks"] = sorted({
                tuple(abs(coordinate) for coordinate in
                      q7.preimage_record(translated, crossing, index)["deck_short"])
                for index, crossing in enumerate(support[1])
            })
            if pairing_ranks and result["mc"]:
                blue_delta = q7.collect_f2(
                    data["precurve"]["delta"] + result["switch"])
                result["ranks"] = {}
                for label, (red_data, red_pre) in earrings.items():
                    certificate = q7.exact_red_blue_homology(
                        red_data, red_pre["vertices"], red_pre["delta"],
                        data, data["precurve"], blue_delta)
                    result["ranks"][label] = certificate["homology"]
        results.append((row, result))
    return results


def f2_span_rank(switches):
    """Rank of finite-support morphisms as F2 coefficient vectors."""
    terms = sorted({term for switch in switches for term in switch}, key=repr)
    positions = {term: index for index, term in enumerate(terms)}
    pivots = {}
    for switch in switches:
        value = sum(1 << positions[term] for term in switch)
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def strict_delta_permutation(first_delta, second_delta, data, vertices):
    """Find an idempotent-preserving labeled-graph permutation, if one exists.

    This is the same restricted notion of strict isomorphism used for the
    explicit S18--S25 relabeling in :mod:`q7_kwz`: a permutation of module
    generators must preserve their skeleton idempotents and carry every
    differential term ``(source, word, target)`` to a term with the identical
    algebra word.  It is deliberately *not* a homotopy-equivalence test.

    Simultaneous directed Weisfeiler--Leman refinement normally makes the
    31-vertex search discrete.  The exact backtracking step handles any color
    classes that remain tied.
    """
    vertices = tuple(vertices)
    if set(first_delta) == set(second_delta):
        return {vertex: vertex for vertex in vertices}

    def graph(delta):
        outgoing = {vertex: [] for vertex in vertices}
        incoming = {vertex: [] for vertex in vertices}
        edge_labels = defaultdict(set)
        for source, word, target in delta:
            if source not in outgoing or target not in incoming:
                raise AssertionError("differential term uses a cancelled vertex")
            outgoing[source].append((word, target))
            incoming[target].append((word, source))
            edge_labels[(source, target)].add(word)
        return outgoing, incoming, edge_labels

    first_out, first_in, first_edges = graph(first_delta)
    second_out, second_in, second_edges = graph(second_delta)
    first_colors = {
        vertex: ("idempotent", data["unique"][vertex]["arc"])
        for vertex in vertices
    }
    second_colors = dict(first_colors)

    while True:
        def signatures(colors, outgoing, incoming):
            return {
                vertex: (
                    colors[vertex],
                    tuple(sorted(
                        (word, colors[target])
                        for word, target in outgoing[vertex])),
                    tuple(sorted(
                        (word, colors[source])
                        for word, source in incoming[vertex])),
                )
                for vertex in vertices
            }

        first_signatures = signatures(first_colors, first_out, first_in)
        second_signatures = signatures(second_colors, second_out, second_in)
        palette = {
            signature: index
            for index, signature in enumerate(sorted(
                set(first_signatures.values()) | set(second_signatures.values()),
                key=repr))
        }
        refined_first = {
            vertex: palette[first_signatures[vertex]] for vertex in vertices
        }
        refined_second = {
            vertex: palette[second_signatures[vertex]] for vertex in vertices
        }
        if Counter(refined_first.values()) != Counter(refined_second.values()):
            return None
        old_classes = len(set(first_colors.values()))
        first_colors, second_colors = refined_first, refined_second
        if len(set(first_colors.values())) == old_classes:
            break

    second_by_color = defaultdict(list)
    for vertex in vertices:
        second_by_color[second_colors[vertex]].append(vertex)
    candidates = {
        vertex: tuple(second_by_color[first_colors[vertex]])
        for vertex in vertices
    }
    mapping = {}
    used = set()

    def compatible(source, image):
        if first_edges.get((source, source), set()) != \
                second_edges.get((image, image), set()):
            return False
        for other, other_image in mapping.items():
            if first_edges.get((source, other), set()) != \
                    second_edges.get((image, other_image), set()):
                return False
            if first_edges.get((other, source), set()) != \
                    second_edges.get((other_image, image), set()):
                return False
        return True

    def search():
        if len(mapping) == len(vertices):
            transported = {
                (mapping[source], word, mapping[target])
                for source, word, target in first_delta
            }
            return mapping.copy() if transported == set(second_delta) else None
        remaining = [vertex for vertex in vertices if vertex not in mapping]
        source = min(
            remaining,
            key=lambda vertex: (
                sum(image not in used for image in candidates[vertex]),
                -(len(first_out[vertex]) + len(first_in[vertex])),
                vertex,
            ))
        for image in candidates[source]:
            if image in used or not compatible(source, image):
                continue
            mapping[source] = image
            used.add(image)
            answer = search()
            if answer is not None:
                return answer
            used.remove(image)
            del mapping[source]
        return None

    return search()


def strict_pair_isomorphism_census(results, blue, selected_pairs):
    """Test the (9,31) two-switch collisions against the S18 object."""
    representatives = {}
    for row, result in results:
        if not result.get("module_matches"):
            continue
        key = frozenset(result["switch"])
        representatives.setdefault(key, (row["name"], result["switch"]))
    by_name = {name: switch for name, switch in representatives.values()}

    translated = q7.translate_to_kwz_chart(blue)
    data = q7.encode_type_d(translated)
    data["precurve"] = q7.geometric_precurve(data)
    vertices = data["precurve"]["vertices"]
    base_delta = data["precurve"]["delta"]
    target_delta = q7.collect_f2(base_delta + by_name["S18"])
    s25_delta = q7.collect_f2(base_delta + by_name["S25"])
    positive_control = strict_delta_permutation(
        target_delta, s25_delta, data, vertices)
    if positive_control is None:
        raise AssertionError(
            "strict-isomorphism search missed the known S18--S25 relabeling")
    known_permutation = {vertex: vertex for vertex in vertices}
    known_permutation.update({4: 22, 22: 4, 5: 21, 21: 5})
    transported = {
        (known_permutation[source], word, known_permutation[target])
        for source, word, target in target_delta
    }
    if transported != set(s25_delta):
        raise AssertionError(
            "known four-vertex S18--S25 relabeling failed the positive control")

    isomorphic = []
    permutations = {}
    for first_name, second_name in selected_pairs:
        candidate_delta = q7.collect_f2(
            base_delta + by_name[first_name] + by_name[second_name])
        permutation = strict_delta_permutation(
            target_delta, candidate_delta, data, vertices)
        if permutation is not None:
            pair = (first_name, second_name)
            isomorphic.append(pair)
            permutations[pair] = permutation
    return isomorphic, permutations


def two_switch_sum_census(results, blue):
    """Compute both closure ranks for every sum of two distinct switches."""
    representatives = {}
    for row, result in results:
        if not result.get("module_matches"):
            continue
        representatives.setdefault(
            frozenset(result["switch"]), (row["name"], result["switch"]))
    if len(representatives) != 38:
        raise AssertionError(
            f"two-switch census expected 38 elements, got {len(representatives)}")

    translated = q7.translate_to_kwz_chart(blue)
    data = q7.encode_type_d(translated)
    data["precurve"] = q7.geometric_precurve(data)
    base_delta = data["precurve"]["delta"]
    earrings = [
        closure.rational_earring_type_d(Fraction(-1, 2)),
        closure.rational_earring_type_d(closure.AUX_SLOPE),
    ]
    rank_census = Counter()
    selected = []
    values = list(representatives.values())
    for (first_name, first), (second_name, second) in itertools.combinations(values, 2):
        switch = q7.collect_f2(first + second)
        residue, _, _ = q7.mc_residue(base_delta, switch)
        if residue:
            raise AssertionError(
                f"{first_name}+{second_name} is not Maurer--Cartan")
        blue_delta = q7.collect_f2(base_delta + switch)
        pair = tuple(
            q7.exact_red_blue_homology(
                red_data, red_pre["vertices"], red_pre["delta"],
                data, data["precurve"], blue_delta)["homology"]
            for red_data, red_pre in earrings
        )
        rank_census[pair] += 1
        if pair == (9, 31):
            selected.append((first_name, second_name))
    return rank_census, selected


def extra_closure_survivors(results, blue, selected_pairs, slopes):
    """Test extra rational-earring ranks on the two-closure collisions."""
    representatives = {}
    for row, result in results:
        if not result.get("module_matches"):
            continue
        key = frozenset(result["switch"])
        representatives.setdefault(key, (row["name"], result["switch"]))
    by_name = {name: switch for name, switch in representatives.values()}
    if "S18" not in by_name:
        raise AssertionError("S18 is absent from the switch representatives")

    translated = q7.translate_to_kwz_chart(blue)
    data = q7.encode_type_d(translated)
    data["precurve"] = q7.geometric_precurve(data)
    base_delta = data["precurve"]["delta"]
    earrings = [closure.rational_earring_type_d(slope) for slope in slopes]

    def ranks(switch):
        blue_delta = q7.collect_f2(base_delta + switch)
        return tuple(
            q7.exact_red_blue_homology(
                red_data, red_pre["vertices"], red_pre["delta"],
                data, data["precurve"], blue_delta)["homology"]
            for red_data, red_pre in earrings
        )

    target = ranks(by_name["S18"])
    survivors = []
    for first_name, second_name in selected_pairs:
        pair_ranks = ranks(q7.collect_f2(
            by_name[first_name] + by_name[second_name]))
        if pair_ranks == target:
            survivors.append((first_name, second_name))
    return target, survivors


def _connector_label(record):
    if record is None:
        return "main"
    seam = "0" if abs(record["seam"]) < 1.0e-8 else "pi"
    return (f"circle={record['circle']} seam={seam} "
            f"{record['lift']}{record['number']} "
            f"angle={math.degrees(record['angle']):.3f}deg")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--switch-census", action="store_true",
        help="encode every connector--main physical smoothing")
    parser.add_argument(
        "--pairing-census", action="store_true",
        help="also compute both exact rational-earring ranks for every switch")
    parser.add_argument(
        "--all-node-census", action="store_true",
        help="test singleton smoothings at all 82 nodes and both closures")
    parser.add_argument(
        "--two-switch-census", action="store_true",
        help="compute both closure ranks for every sum of two distinct switches")
    parser.add_argument(
        "--strict-pair-census", action="store_true",
        help="with --two-switch-census, test (9,31) sums for strict relabeling")
    parser.add_argument(
        "--extra-slope", action="append", type=Fraction, default=[],
        help="with --two-switch-census, test an additional rational closure")
    args = parser.parse_args(argv)

    _, blue, xinfo = q7.build_q7()
    rows = provenance_census(blue, xinfo)
    counts = Counter(row["kind"] for row in rows)
    print("== q=7 C3 provenance census ==")
    print(f"crossing orbits={len(rows)}; provenance={dict(sorted(counts.items()))}")
    if len(rows) != 82 or dict(counts) != EXPECTED_PROVENANCE:
        raise AssertionError("q=7 crossing/provenance census is not stable")

    canonical = canonical_connector_orbits(rows, xinfo)
    if canonical != EXPECTED_CANONICAL:
        raise AssertionError(
            f"canonical connector nodes changed: {canonical} != {EXPECTED_CANONICAL}")
    print("canonical connector--connector nodes by fibre circle:")
    for circle, name in sorted(canonical.items()):
        print(f"  circle {circle}: {name}")

    named = locate_named(rows)
    if any(named[name]["name"] != name for name in q7.TARGETS):
        raise AssertionError("a named q=7 target moved to a different crossing orbit")
    print("named physical switches:")
    for name in ("S18", "S25", "S69", "S74"):
        row = named[name]
        print(f"  {name} -> {row['name']} at "
              f"({row['point'][0]:.7f},{row['point'][1]:.7f}) "
              f"[{row['kind']}]")
        for lift in row["lifts"]:
            print(f"    edges={lift['edges']}: "
                  f"A={_connector_label(lift['A'])}; "
                  f"B={_connector_label(lift['B'])}")
        if row["kind"] != "connector--main":
            raise AssertionError(f"{name} is not connector--main")
    if set(canonical.values()) & {row["name"] for row in named.values()}:
        raise AssertionError("a named switch was mistaken for a canonical connector node")

    if args.extra_slope and not args.two_switch_census:
        parser.error("--extra-slope requires --two-switch-census")
    if args.strict_pair_census and not args.two_switch_census:
        parser.error("--strict-pair-census requires --two-switch-census")
    if args.all_node_census and args.two_switch_census:
        parser.error("--all-node-census is a singleton census; run the "
                     "connector--main two-switch census separately")

    run_switches = (args.switch_census or args.pairing_census
                    or args.all_node_census or args.two_switch_census)
    if run_switches:
        scope = "all-node" if args.all_node_census else "connector--main"
        print(f"\n== module-preserving {scope} physical-switch census ==")
        pairing_ranks = (args.pairing_census or args.all_node_census
                         or args.two_switch_census)
        results = switch_census(
            rows, blue, pairing_ranks=pairing_ranks,
            all_nodes=args.all_node_census)
        module_preserving = []
        mc_switches = []
        closed_switches = []
        square_zero_switches = []
        delta_homogeneous = []
        selected_by_ranks = []
        excluded = []
        switch_fibres = defaultdict(list)
        rank_pairs = {}
        for row, result in results:
            if not result["module_matches"]:
                excluded.append(row["name"])
                print(f"  {row['name']}: excluded ({result['error']})")
                continue
            module_preserving.append(row["name"])
            # A morphism is an F2 set of terms.  ``collect_f2`` preserves
            # insertion order, so tuple keys would incorrectly distinguish
            # the S0 and S14 presentations of the same four-arrow element.
            switch_fibres[frozenset(result["switch"])].append(row["name"])
            if result["mc"]:
                mc_switches.append(row["name"])
            if result["d_terms"] == 0:
                closed_switches.append(row["name"])
            if result["square_terms"] == 0:
                square_zero_switches.append(row["name"])
            if result["delta_homogeneous"]:
                delta_homogeneous.append(row["name"])
            quantum = sorted({str(degree[0]) for degree in result["degrees"]})
            rank_text = ""
            if pairing_ranks and result["mc"]:
                pair = (result["ranks"]["original"],
                        result["ranks"]["auxiliary"])
                rank_pairs[row["name"]] = pair
                rank_text = f" ranks={pair}"
                if pair == (9, 31):
                    selected_by_ranks.append(row["name"])
            print(f"  {row['name']}: terms={result['terms']} "
                  f"MC={result['mc']} delta-hom={result['delta_homogeneous']} "
                  f"q={quantum} residue={result['residue_terms']}" + rank_text)
        duplicate_switches = {
            frozenset(names) for names in switch_fibres.values() if len(names) > 1
        }
        switch_span_rank = f2_span_rank(
            [list(switch) for switch in switch_fibres])
        print(f"{scope} orbits tested={len(results)}")
        print(f"distinct switch elements={len(switch_fibres)}")
        print(f"F2 switch-span rank={switch_span_rank}")
        print(f"module-preserving={module_preserving}")
        print(f"Maurer--Cartan switches={mc_switches}")
        print(f"D-closed switches={closed_switches}")
        print(f"square-zero switches={square_zero_switches}")
        print(f"delta-degree -1 switches={delta_homogeneous}")
        if pairing_ranks:
            node_rank_census = Counter(rank_pairs.values())
            distinct_rank_census = Counter()
            for names in switch_fibres.values():
                ranks = {rank_pairs[name] for name in names}
                if len(ranks) != 1:
                    raise AssertionError(
                        f"duplicate switch has inconsistent ranks: {names}")
                distinct_rank_census.update(ranks)
            print(f"node rank-pair census={dict(sorted(node_rank_census.items()))}")
            print("distinct-switch rank-pair census="
                  f"{dict(sorted(distinct_rank_census.items()))}")
            print(f"rank-pair (9,31) switches={selected_by_ranks}")

        expected_count = 82 if args.all_node_census else 52
        expected_exclusions = (EXPECTED_ALL_MODULE_EXCLUSIONS
                               if args.all_node_census
                               else EXPECTED_MODULE_EXCLUSIONS)
        expected_admissible = 73 if args.all_node_census else 46
        expected_distinct = 45 if args.all_node_census else 38
        expected_duplicates = (EXPECTED_ALL_DUPLICATE_SWITCHES
                               if args.all_node_census
                               else EXPECTED_DUPLICATE_SWITCHES)
        expected_span_rank = 45 if args.all_node_census else 38
        expected_node_ranks = (EXPECTED_ALL_NODE_RANK_CENSUS
                               if args.all_node_census
                               else EXPECTED_NODE_RANK_CENSUS)
        expected_distinct_ranks = (EXPECTED_ALL_DISTINCT_RANK_CENSUS
                                   if args.all_node_census
                                   else EXPECTED_DISTINCT_RANK_CENSUS)
        if len(results) != expected_count:
            raise AssertionError(
                f"{scope} census no longer has {expected_count} orbits")
        if set(excluded) != expected_exclusions:
            raise AssertionError(
                f"module exclusions changed: {excluded} != "
                f"{sorted(expected_exclusions)}")
        if len(module_preserving) != expected_admissible:
            raise AssertionError(
                f"physical-switch census no longer has {expected_admissible} "
                "admissible nodes")
        if mc_switches != module_preserving:
            raise AssertionError("a module-preserving physical switch failed MC")
        if closed_switches != module_preserving:
            raise AssertionError("a module-preserving physical switch was not D-closed")
        if square_zero_switches != module_preserving:
            raise AssertionError("a module-preserving physical switch had nonzero square")
        if delta_homogeneous != module_preserving:
            raise AssertionError("a module-preserving physical switch changed delta degree")
        if (len(switch_fibres) != expected_distinct
                or duplicate_switches != expected_duplicates):
            raise AssertionError(
                f"the {expected_admissible} physical nodes no longer give "
                f"the expected {expected_distinct} switch elements")
        if switch_span_rank != expected_span_rank:
            raise AssertionError(
                f"the {expected_distinct} distinct physical switches are "
                "not independent")
        if pairing_ranks:
            if node_rank_census != expected_node_ranks:
                raise AssertionError(
                    "the node-orbit two-closure rank census changed")
            if distinct_rank_census != expected_distinct_ranks:
                raise AssertionError(
                    "the distinct-switch two-closure rank census changed")
            if selected_by_ranks != ["S18", "S25"]:
                raise AssertionError("the (9,31) closure selector is no longer unique")
            if args.all_node_census:
                deck_selected = [
                    row["name"] for row, result in results
                    if result.get("module_matches")
                    and (2, 2) in result["short_decks"]
                ]
                print(f"short-deck (2,2) switches={deck_selected}")
                if deck_selected != ["S18", "S25"]:
                    raise AssertionError(
                        "the all-node short-deck (2,2) selector changed")
            expected_named_ranks = {
                "S18": (9, 31),
                "S25": (9, 31),
                "S69": (9, 23),
                "S74": (9, 25),
            }
            if {name: rank_pairs[name] for name in expected_named_ranks} != expected_named_ranks:
                raise AssertionError("the four named switch rank pairs changed")

        if args.two_switch_census:
            print("\n== all sums of two distinct switch elements ==")
            two_rank_census, two_selected = two_switch_sum_census(results, blue)
            print(f"rank-pair census={dict(sorted(two_rank_census.items()))}")
            print(f"rank-pair (9,31) two-switch sums={two_selected}")
            if sum(two_rank_census.values()) != 703:
                raise AssertionError("two-switch census did not test C(38,2)=703 sums")
            if two_rank_census != EXPECTED_TWO_SWITCH_RANK_CENSUS:
                raise AssertionError("the two-switch two-closure rank census changed")
            if len(two_selected) != 41:
                raise AssertionError("the two-switch (9,31) multiplicity changed")
            if args.strict_pair_census:
                isomorphic, permutations = strict_pair_isomorphism_census(
                    results, blue, two_selected)
                print("strictly permutation-isomorphic to S18="
                      f"{isomorphic}")
                print("nontrivial permutation supports=" + str({
                    pair: sorted(
                        vertex for vertex, image in permutation.items()
                        if vertex != image)
                    for pair, permutation in permutations.items()
                }))
                if isomorphic:
                    raise AssertionError(
                        "a two-switch (9,31) collision became a relabeling of S18")
            if args.extra_slope:
                target, survivors = extra_closure_survivors(
                    results, blue, two_selected, args.extra_slope)
                print(f"extra slopes={args.extra_slope}; S18 ranks={target}")
                print(f"two-switch survivors after extra slopes={survivors}")

    print("\nPROVED BY THIS CERTIFICATE:")
    if run_switches:
        if args.all_node_census:
            print("  provenance and singleton reconnection data for all 82 crossings")
        else:
            print("  provenance of all 82 crossings and the exhaustive connector--"
                  "main physical-reconnection census")
    else:
        print("  provenance of all 82 crossings")
    print("NOT PROVED (Gate-1 analytic localization lemma):")
    print("  every rigid C3 figure-eight output is represented by one "
          "module-preserving physical reconnection or by the named span")
    return 0


if __name__ == "__main__":
    sys.exit(main())
