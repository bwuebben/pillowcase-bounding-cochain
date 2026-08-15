#!/usr/bin/env python3
"""Compare the corrected analytic C3 q=7 curve with the finite PL model.

This is a numerical continuation and exact finite-algebra regression for Gate 1.
It composes Smith's corrected C3 correspondence with Q_{1/3} and Q_{1/7},
traces one (21,10) lift, and checks three coordinate-independent records:

* the short deck class at every transverse double point;
* the signed face-word interval of the admissible KWZ two-arc encoding.
* optionally, the exact two-closure ranks of every module-preserving singleton
  smoothing and the mixed Maurer--Cartan products between those switches.

The floating-point trace is not an interval-arithmetic proof.  Nor does this
script prove the CHKK instanton--pillowcase assignment.  It does show that the
continued group-word model and ``resolve.py`` encode strictly isomorphic
31-generator type-D base objects, and that the two S18/S25 deck-(2,2) nodes are
present on the continued curve rather than artifacts of the sinusoidal PL
connectors.

Run from the repository root or from ``pillowcase``:

    python3 pillowcase/c3_q7_compare.py
    python3 pillowcase/c3_q7_compare.py --stability
    python3 pillowcase/c3_q7_compare.py --stability --switch-census
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import c3_perturbed as c3
import q7_kwz as q7
import q7_closure_probe as closure
import q7_quilt_census as quilt
from bigons import _tdist
from earring import P_point
from maurer_cartan import orbit_group
from polygons import self_intersections_detailed
from surgery_check import smooth_curve


EXPECTED_DECK_CENSUS = Counter({
    (3, 1): 8,
    (10, 5): 8,
    (8, 4): 6,
    (5, 3): 4,
    (6, 2): 4,
    (5, 2): 2,
    (6, 3): 2,
    (4, 2): 2,
    (7, 3): 2,
    (3, 2): 2,
    (2, 1): 2,
    (8, 3): 2,
    (9, 4): 2,
    (2, 2): 2,
    (1, 1): 2,
})

EXPECTED_SINGLETON_RANK_CENSUS = Counter({
    (5, 23): 3,
    (5, 25): 2,
    (7, 23): 6,
    (7, 25): 21,
    (7, 27): 2,
    (9, 23): 1,
    (9, 25): 6,
    (9, 31): 2,
})


def canonical_deck(vector):
    """Choose the sign whose first nonzero coordinate is positive."""
    first, second = vector
    if first < 0 or (first == 0 and second < 0):
        return (-first, -second)
    return vector


def interpolate_state(lift, edge, parameter):
    return tuple(
        lift[edge][index]
        + parameter * (lift[edge + 1][index] - lift[edge][index])
        for index in range(4)
    )


def normalized_crossing_angle(crossing):
    first, second = crossing["dirA"], crossing["dirB"]
    denominator = math.hypot(*first) * math.hypot(*second)
    if denominator == 0.0:
        return 0.0
    return abs(first[0] * second[1] - first[1] * second[0]) / denominator


def iota_crossing_pairs(crossings):
    """Pair detailed T2 crossings by iota, with a nearest-neighbor margin gate."""
    partners = {}
    nearest_distances = {}
    margins = {}
    for index, crossing in enumerate(crossings):
        target = ((-crossing["pt"][0]) % c3.TAU,
                  (-crossing["pt"][1]) % c3.TAU)
        ranked = sorted(
            (_tdist(target, candidate["pt"]), other)
            for other, candidate in enumerate(crossings) if other != index
        )
        if len(ranked) < 2:
            raise AssertionError("not enough crossings to certify an iota pair")
        nearest, partner = ranked[0]
        runner_up = ranked[1][0]
        partners[index] = partner
        nearest_distances[index] = nearest
        margins[index] = runner_up / max(nearest, 1.0e-15)

    if any(partners.get(partner) != index
           for index, partner in partners.items()):
        raise AssertionError("nearest iota partners are not mutual")
    if max(nearest_distances.values()) >= 1.0e-3:
        raise AssertionError("polyline is not sufficiently iota-symmetric")
    if min(margins.values()) <= 8.0:
        raise AssertionError("an iota partner is not separated from other crossings")

    pairs = [(index, partner) for index, partner in partners.items()
             if index < partner]
    return pairs, {
        "worst_match": max(nearest_distances.values()),
        "minimum_margin": min(margins.values()),
    }


def analytic_crossing_records(curve, lift):
    crossings = self_intersections_detailed(curve)
    if len(crossings) != 100:
        raise AssertionError(
            f"corrected q=7 curve has {len(crossings)} rather than 100 T2 crossings")
    pairs, pairing_diagnostics = iota_crossing_pairs(crossings)
    if len(pairs) != 50:
        raise AssertionError("corrected q=7 crossings do not form 50 iota pairs")

    total_deck = (21, 10)
    maximum_deck_error = 0.0
    records = []
    for first_index, second_index in pairs:
        crossing = crossings[first_index]
        first = interpolate_state(lift, crossing["kA"], crossing["tA"])
        second = interpolate_state(lift, crossing["kB"], crossing["tB"])
        raw = ((second[0] - first[0]) / c3.TAU,
               (second[1] - first[1]) / c3.TAU)
        forward = tuple(round(value) for value in raw)
        maximum_deck_error = max(
            maximum_deck_error,
            max(abs(value - rounded) for value, rounded in zip(raw, forward)),
        )
        reverse = (total_deck[0] - forward[0],
                   total_deck[1] - forward[1])
        short = min((forward, reverse),
                    key=lambda value: (abs(value[0]) + abs(value[1]), value))
        records.append({
            "point": P_point(crossing["pt"]),
            "deck": canonical_deck(short),
            "angle_sine": min(
                normalized_crossing_angle(crossings[first_index]),
                normalized_crossing_angle(crossings[second_index]),
            ),
        })

    if maximum_deck_error >= 2.0e-3:
        raise AssertionError("a crossing deck difference is not near an integer")
    census = Counter(record["deck"] for record in records)
    if census != EXPECTED_DECK_CENSUS:
        raise AssertionError(
            f"analytic deck census changed: {census} != {EXPECTED_DECK_CENSUS}")
    minimum_angle_sine = min(record["angle_sine"] for record in records)
    if minimum_angle_sine <= 1.0e-4:
        raise AssertionError("analytic trace contains a numerically near-tangent crossing")
    pairing_diagnostics.update({
        "maximum_deck_error": maximum_deck_error,
        "minimum_angle_sine": minimum_angle_sine,
    })
    return records, pairing_diagnostics


def coalesce_iota_skeleton_hits(data):
    """Impose the exact double-cover quotient on approximate skeleton heights.

    Near the stereographic point at infinity, a tiny T2 symmetry error is
    amplified in the affine height.  The two occurrences of each geometric
    dot are nevertheless consecutive by height and separated from neighboring
    dots by a large relative gap.  This routine checks that separation before
    replacing each pair by its mean height.
    """
    maximum_gap = 0.0
    maximum_relative_gap = 0.0
    for arc in ("L", "R"):
        indices = sorted(
            (index for index, hit in enumerate(data["hits"])
             if hit["arc"] == arc),
            key=lambda index: data["hits"][index]["y"],
        )
        if len(indices) % 2:
            raise AssertionError(f"{arc} skeleton has an odd number of T2 hits")
        for offset in range(0, len(indices), 2):
            first, second = indices[offset:offset + 2]
            y_first = data["hits"][first]["y"]
            y_second = data["hits"][second]["y"]
            pair_gap = y_second - y_first
            neighbor_gaps = []
            if offset > 0:
                neighbor_gaps.append(
                    y_first - data["hits"][indices[offset - 1]]["y"])
            if offset + 2 < len(indices):
                neighbor_gaps.append(
                    data["hits"][indices[offset + 2]]["y"] - y_second)
            neighbor_gap = min(neighbor_gaps, default=math.inf)
            relative = pair_gap / neighbor_gap if math.isfinite(neighbor_gap) else 0.0
            if relative >= 0.05:
                raise AssertionError(
                    f"{arc} iota-hit pair is not isolated: ratio={relative:.3e}")
            mean = 0.5 * (y_first + y_second)
            data["hits"][first]["y"] = mean
            data["hits"][second]["y"] = mean
            maximum_gap = max(maximum_gap, pair_gap)
            maximum_relative_gap = max(maximum_relative_gap, relative)

    data["unique"] = q7.unique_hits(data["hits"], ytol=1.0e-9)
    return {
        "maximum_pair_gap": maximum_gap,
        "maximum_relative_gap": maximum_relative_gap,
    }


def average_paired_component_hits(first, second):
    """Identify the skeleton hits of two components exchanged by iota.

    The numerical trace is only approximately iota-invariant.  The two
    components therefore have slightly different affine skeleton heights.
    Corresponding hits have the same arc and order; the next geometric dot
    gives a scale-free separation gate before the two heights are averaged.
    """
    maximum_gap = 0.0
    maximum_relative_gap = 0.0
    for arc in ("L", "R"):
        first_indices = sorted(
            (index for index, hit in enumerate(first["hits"])
             if hit["arc"] == arc),
            key=lambda index: first["hits"][index]["y"],
        )
        second_indices = sorted(
            (index for index, hit in enumerate(second["hits"])
             if hit["arc"] == arc),
            key=lambda index: second["hits"][index]["y"],
        )
        if len(first_indices) != len(second_indices):
            raise AssertionError(
                f"iota-paired {arc} components have different hit counts")
        means = [
            0.5 * (first["hits"][i]["y"] + second["hits"][j]["y"])
            for i, j in zip(first_indices, second_indices)
        ]
        for rank, (i, j) in enumerate(zip(first_indices, second_indices)):
            first_y = first["hits"][i]["y"]
            second_y = second["hits"][j]["y"]
            gap = abs(first_y - second_y)
            neighbors = []
            if rank:
                neighbors.append(means[rank] - means[rank - 1])
            if rank + 1 < len(means):
                neighbors.append(means[rank + 1] - means[rank])
            scale = min(neighbors, default=math.inf)
            relative = gap / scale if math.isfinite(scale) else 0.0
            if relative >= 0.05:
                raise AssertionError(
                    f"{arc} component pair is not isolated: ratio={relative:.3e}")
            first["hits"][i]["y"] = means[rank]
            maximum_gap = max(maximum_gap, gap)
            maximum_relative_gap = max(maximum_relative_gap, relative)
    first["unique"] = q7.unique_hits(first["hits"], ytol=1.0e-9)
    if len(first["unique"]) != len(first["hits"]):
        raise AssertionError("one component has repeated quotient skeleton dots")
    return {
        "maximum_pair_gap": maximum_gap,
        "maximum_relative_gap": maximum_relative_gap,
    }


def match_component_dots_by_order(records, base_data):
    """Match a quotient collection to the base dots by arc and height order."""
    matching = {}
    maximum_relative_displacement = 0.0
    for arc in ("L", "R"):
        local = sorted(
            (data["unique"][index]["y"], record, index)
            for record, (data, _) in enumerate(records)
            for index in range(len(data["unique"]))
            if data["unique"][index]["arc"] == arc
        )
        base = sorted(
            (dot["y"], index)
            for index, dot in enumerate(base_data["unique"])
            if dot["arc"] == arc
        )
        if len(local) != len(base):
            raise AssertionError(
                f"smoothed {arc} dot count {len(local)} != base {len(base)}")
        for rank, ((height, record, index), (base_height, base_index)) in enumerate(
                zip(local, base)):
            neighbors = []
            if rank:
                neighbors.append(base_height - base[rank - 1][0])
            if rank + 1 < len(base):
                neighbors.append(base[rank + 1][0] - base_height)
            scale = min(neighbors, default=math.inf)
            relative = (abs(height - base_height) / scale
                        if math.isfinite(scale) else 0.0)
            if relative >= 0.05:
                raise AssertionError(
                    f"smoothed {arc} dot changed order: ratio={relative:.3e}")
            maximum_relative_displacement = max(
                maximum_relative_displacement, relative)
            matching[(record, index)] = base_index
    return matching, maximum_relative_displacement


def directional_iota_pairing(first, second, source_pairing=1):
    """Transport a smoothing sector through d(iota)=-I using half-rays."""
    chord_sets = {
        1: (("A-", "B+"), ("B-", "A+")),
        2: (("A-", "B-"), ("A+", "B+")),
    }

    def unit(vector):
        norm = math.hypot(*vector)
        if norm == 0.0:
            raise AssertionError("crossing branch has zero tangent")
        return vector[0] / norm, vector[1] / norm

    def rays(crossing):
        branch_a = unit(crossing["dirA"])
        branch_b = unit(crossing["dirB"])
        return {
            "A-": (-branch_a[0], -branch_a[1]),
            "A+": branch_a,
            "B-": (-branch_b[0], -branch_b[1]),
            "B+": branch_b,
        }

    first_rays = rays(first)
    second_rays = rays(second)
    labels = tuple(first_rays)
    costs = {1: math.inf, 2: math.inf}
    target_matchings = {
        pairing: {frozenset(chord) for chord in chords}
        for pairing, chords in chord_sets.items()
    }
    for permutation in itertools.permutations(labels):
        matching = dict(zip(labels, permutation))
        image = {
            frozenset((matching[start], matching[end]))
            for start, end in chord_sets[source_pairing]
        }
        target_pairing = next(
            (pairing for pairing, target in target_matchings.items()
             if image == target), None)
        if target_pairing is None:
            continue
        cost = sum(
            (-first_rays[label][0]
             - second_rays[matching[label]][0]) ** 2
            + (-first_rays[label][1]
               - second_rays[matching[label]][1]) ** 2
            for label in labels
        )
        costs[target_pairing] = min(costs[target_pairing], cost)
    pairing = min(costs, key=costs.get)
    alternative = 3 - pairing
    margin = costs[alternative] / max(costs[pairing], 1.0e-15)
    if margin <= 100.0:
        raise AssertionError(
            f"iota smoothing sector is ambiguous: costs={costs}")
    return pairing, math.sqrt(costs[pairing] / 4.0), margin


def analytic_physical_switch(blue, support, base_data, eps=0.002,
                             subdivisions=40):
    """Encode one actual-C3 smoothing on the fixed 31-generator module."""
    point, lifts = support
    if len(lifts) != 2:
        return {"module_matches": False,
                "error": f"{len(lifts)} torus lifts"}
    pairing, direction_error, pairing_margin = directional_iota_pairing(
        lifts[0], lifts[1])
    components = smooth_curve(
        blue, [(lifts[0], 1), (lifts[1], pairing)], eps)
    component_data = [
        q7.encode_type_d(component, n_samp=subdivisions)
        for component in components
    ]
    quotient_diagnostics = []
    selected = []
    if len(components) == 1:
        quotient_diagnostics.append(
            coalesce_iota_skeleton_hits(component_data[0]))
        selected.append(component_data[0])
    elif len(components) == 3:
        decks = [q7.deck_of_path(component) for component in components]
        paired = [
            (first, second)
            for first in range(3) for second in range(first + 1, 3)
            if decks[first] == decks[second]
        ]
        if len(paired) != 1:
            return {"module_matches": False,
                    "error": f"component decks do not give one iota pair: {decks}"}
        first, second = paired[0]
        invariant = next(index for index in range(3)
                         if index not in (first, second))
        quotient_diagnostics.append(
            coalesce_iota_skeleton_hits(component_data[invariant]))
        quotient_diagnostics.append(
            average_paired_component_hits(
                component_data[first], component_data[second]))
        selected.extend((component_data[invariant], component_data[first]))
    else:
        return {"module_matches": False,
                "error": f"smoothing produced {len(components)} components"}

    records = []
    for data in selected:
        invariant = len(data["hits"]) == 2 * len(data["unique"])
        precurve = (q7.geometric_precurve(data, ytol=1.0e-9)
                    if invariant else q7.cyclic_precurve(data))
        records.append((data, precurve))
    try:
        matching, displacement = match_component_dots_by_order(
            records, base_data)
    except AssertionError as error:
        return {"module_matches": False, "error": str(error)}

    vertices = []
    smoothed_delta = []
    for record, (_, precurve) in enumerate(records):
        vertices.extend(
            matching[(record, vertex)] for vertex in precurve["vertices"])
        smoothed_delta.extend(
            (matching[(record, source)], word, matching[(record, target)])
            for source, word, target in precurve["delta"]
        )
    base = base_data["precurve"]
    if (len(vertices) != len(set(vertices))
            or set(vertices) != set(base["vertices"])):
        return {
            "module_matches": False,
            "error": ("cancelled module changed: "
                      f"{len(vertices)} vertices, {len(set(vertices))} distinct, "
                      f"{len(set(vertices) & set(base['vertices']))} base"),
        }
    smoothed_delta = q7.collect_f2(smoothed_delta)
    switch = q7.collect_f2(base["delta"] + smoothed_delta)
    residue, differential, square = q7.mc_residue(base["delta"], switch)
    return {
        "module_matches": True,
        "switch": switch,
        "terms": len(switch),
        "mc": not residue,
        "d_terms": len(differential),
        "square_terms": len(square),
        "residue_terms": len(residue),
        "iota_direction_error": direction_error,
        "iota_pairing_margin": pairing_margin,
        "maximum_relative_dot_displacement": displacement,
        "maximum_relative_quotient_gap": max(
            item["maximum_relative_gap"] for item in quotient_diagnostics),
        "component_decks": [
            q7.deck_of_path(component) for component in components
        ],
        "point": point,
    }


def analytic_precurve(curve, subdivisions):
    translated = q7.translate_to_kwz_chart(curve)
    data = q7.encode_type_d(translated, n_samp=subdivisions)
    quotient_diagnostics = coalesce_iota_skeleton_hits(data)
    precurve = q7.geometric_precurve(data, ytol=1.0e-9)
    data["precurve"] = precurve
    expected = (
        len(data["hits"]), len(data["unique"]), len(precurve["path"]),
        len(precurve["cancellations"]), len(precurve["vertices"]),
        len(precurve["delta"]), len(precurve["residue"]),
    )
    if expected != (74, 37, 37, 3, 31, 30, 0):
        raise AssertionError(
            "analytic KWZ reduction changed: "
            f"{expected} != (74,37,37,3,31,30,0)")
    if not precurve["paired_ok"]:
        raise AssertionError("analytic geometric joins do not pair on the T2 cover")
    return data, precurve, quotient_diagnostics


def base_precurve():
    _, blue, xinfo = q7.build_q7()
    data = q7.encode_type_d(q7.translate_to_kwz_chart(blue))
    precurve = q7.geometric_precurve(data)
    if (len(precurve["vertices"]), len(precurve["delta"]),
            len(precurve["residue"])) != (31, 30, 0):
        raise AssertionError("PL base object is no longer the 31/30 type-D module")
    return blue, xinfo, data, precurve


def strict_interval_isomorphism(analytic, base):
    analytic_words = tuple(analytic["words"])
    expected_words = tuple(
        (face, -length) for face, length in reversed(base["words"]))
    if analytic_words != expected_words:
        raise AssertionError("analytic and PL signed face-word intervals differ")
    matching = {
        analytic["path"][index]: base["path"][-1 - index]
        for index in range(len(analytic["path"]))
    }
    mapped_delta = q7.collect_f2(
        (matching[source], word, matching[target])
        for source, word, target in analytic["delta"]
    )
    if set(mapped_delta) != set(base["delta"]):
        raise AssertionError("the interval relabeling does not intertwine delta")
    return matching


def pl_deck_two_nodes(blue, xinfo):
    rows = quilt.provenance_census(blue, xinfo)
    names = []
    for row in rows:
        if row["kind"] != "connector--main":
            continue
        record = q7.preimage_record(blue, row["preimages"][0], 0)
        if canonical_deck(record["deck_short"]) == (2, 2):
            names.append(row["name"])
    if names != ["S18", "S25"]:
        raise AssertionError(f"PL deck-(2,2) nodes changed: {names}")
    return names


def _switch_key(switch):
    return frozenset(switch)


def _f2_span_rank(switches):
    """Rank of finite-support morphisms viewed as F2 coefficient vectors."""
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


def analytic_switch_census(trace_result, base):
    """Census singleton surgeries at every crossing of the actual C3 trace."""
    pl_blue, pl_xinfo, pl_data, pl_precurve = base
    analytic_curve = q7.translate_to_kwz_chart(trace_result["curve"])
    supports = orbit_group(analytic_curve)
    if len(supports) != 50 or any(len(lifts) != 2 for _, lifts in supports):
        raise AssertionError("analytic switch census is not 50 two-lift orbits")

    pl_lookup = defaultdict(list)
    for row, result in quilt.switch_census(
            quilt.provenance_census(pl_blue, pl_xinfo), pl_blue):
        if result.get("module_matches"):
            pl_lookup[_switch_key(result["switch"])].append(row["name"])

    earrings = {
        "original": closure.rational_earring_type_d(Fraction(-1, 2)),
        "auxiliary": closure.rational_earring_type_d(closure.AUX_SLOPE),
    }
    rank_cache = {}

    def rank_pair(mapped_switch):
        key = _switch_key(mapped_switch)
        if key not in rank_cache:
            blue_delta = q7.collect_f2(pl_precurve["delta"] + mapped_switch)
            rank_cache[key] = tuple(
                q7.exact_red_blue_homology(
                    red_data, red_precurve["vertices"],
                    red_precurve["delta"], pl_data, pl_precurve,
                    blue_delta)["homology"]
                for red_data, red_precurve in earrings.values()
            )
        return rank_cache[key]

    rows = []
    for index, support in enumerate(supports):
        original_point = (
            (support[0][0] + q7.KWZ_SPECIAL[0]) % c3.TAU,
            (support[0][1] + q7.KWZ_SPECIAL[1]) % c3.TAU,
        )
        deck = canonical_deck(
            q7.preimage_record(
                analytic_curve, support[1][0], index)["deck_short"])
        try:
            result = analytic_physical_switch(
                analytic_curve, support, trace_result["switch_data"])
        except (AssertionError, KeyError, ValueError) as error:
            result = {"module_matches": False,
                      "error": f"{type(error).__name__}: {error}"}
        row = {
            "index": index,
            "point": P_point(original_point),
            "deck": deck,
            "result": result,
        }
        if result.get("module_matches"):
            if result["terms"] != 4 or not result["mc"]:
                raise AssertionError(
                    f"analytic node A{index} is not a four-arrow MC switch")
            mapped = q7.collect_f2(
                (trace_result["switch_matching"][source], word,
                 trace_result["switch_matching"][target])
                for source, word, target in result["switch"]
            )
            row["mapped_switch"] = mapped
            row["pl_names"] = tuple(pl_lookup.get(_switch_key(mapped), ()))
            row["ranks"] = rank_pair(mapped)
        rows.append(row)

    module_rows = [row for row in rows
                   if row["result"].get("module_matches")]
    excluded = [row for row in rows
                if not row["result"].get("module_matches")]
    distinct = {_switch_key(row["mapped_switch"]) for row in module_rows}
    rank_census = Counter(row["ranks"] for row in module_rows)
    selected = [row for row in module_rows if row["ranks"] == (9, 31)]
    deck_two = [row for row in module_rows if row["deck"] == (2, 2)]
    switches = [row["mapped_switch"] for row in module_rows]
    mixed_residues = []
    for first_index, first in enumerate(switches):
        for second_index, second in enumerate(
                switches[first_index + 1:], start=first_index + 1):
            residue = q7.collect_f2(
                q7.apply_mul(first, second) + q7.apply_mul(second, first))
            if residue:
                mixed_residues.append((first_index, second_index, residue))
    switch_span_rank = _f2_span_rank(switches)

    print("\n== actual-C3 singleton switch census ==")
    print(f"analytic orbits=50; module-preserving={len(module_rows)}; "
          f"presentation-changing={len(excluded)}; "
          f"distinct switches={len(distinct)}")
    print("two-closure rank-pair census: "
          + ", ".join(
              f"{pair}:{count}" for pair, count in sorted(rank_census.items())))
    print("mixed MC anticommutators=0; "
          f"all 2^{switch_span_rank} distinct F2 switch sums are MC")
    print("analytic switches absent from the PL connector census="
          f"{sum(not row['pl_names'] for row in module_rows)}")
    print("rank-pair (9,31) analytic nodes:")
    for row in selected:
        print(f"  A{row['index']}: deck={row['deck']}; "
              f"PL={','.join(row['pl_names']) or 'none'}")
    print("presentation-changing singleton surgeries:")
    for row in excluded:
        print(f"  A{row['index']}: deck={row['deck']}; "
              f"{row['result']['error']}")

    if len(module_rows) != 43 or len(excluded) != 7 or len(distinct) != 43:
        raise AssertionError("actual-C3 module-preserving census changed")
    if rank_census != EXPECTED_SINGLETON_RANK_CENSUS:
        raise AssertionError(
            "actual-C3 singleton two-closure rank census changed")
    if mixed_residues:
        raise AssertionError(
            "two actual-C3 switches have a nonzero mixed MC residue")
    if switch_span_rank != 43:
        raise AssertionError("actual-C3 switches are no longer independent")
    if ({row["deck"] for row in selected} != {(2, 2)}
            or {name for row in selected for name in row["pl_names"]}
               != {"S18", "S25"}
            or len(selected) != 2
            or len(deck_two) != 2):
        raise AssertionError(
            "the analytic rank-pair selector no longer isolates S18/S25")

    return {
        "rows": rows,
        "rank_census": rank_census,
        "selected": selected,
        "switch_span_rank": switch_span_rank,
        "signature": tuple(sorted(
            (row["deck"], tuple(sorted(row["mapped_switch"], key=repr)),
             row["ranks"])
            for row in module_rows
        )),
    }


def run_one(perturbation, step, base):
    blue, _, _, base_pre = base
    curve, lift, trace = c3.trace_corrected_main_arclength(
        perturbation, step=step, max_steps=12000)
    if trace["maximum_residual"] > 2.1e-10:
        raise AssertionError("pseudo-arclength residual exceeded its gate")
    if trace["maximum_tangent_residual"] > 3.0e-7:
        raise AssertionError("computed tangent does not annihilate every equation")

    crossing_records, crossing_diagnostics = analytic_crossing_records(curve, lift)
    deck_two = sorted(record["point"] for record in crossing_records
                      if record["deck"] == (2, 2))
    if len(deck_two) != 2 or deck_two[0][0] >= c3.PI / 2.0 \
            or deck_two[1][0] <= c3.PI / 2.0:
        raise AssertionError("the two analytic deck-(2,2) nodes changed sides")

    signatures = []
    quotient_diagnostics = []
    matching = None
    switch_data = None
    switch_matching = None
    for subdivisions in (20, 40, 80):
        data, precurve, quotient = analytic_precurve(curve, subdivisions)
        matching = strict_interval_isomorphism(precurve, base_pre)
        if subdivisions == 40:
            switch_data = data
            switch_matching = matching
        signatures.append((tuple(precurve["words"]), tuple(precurve["delta"])))
        quotient_diagnostics.append(quotient)
    if not all(signature == signatures[0] for signature in signatures[1:]):
        raise AssertionError("analytic KWZ encoding depends on root subdivisions")

    print(f"\n== corrected C3 q=7 trace: t={perturbation:g} ==")
    print(f"points={len(curve) - 1}; steps={trace['steps']}; "
          f"raw closure error={trace['raw_closure_error']:.3e}")
    print(f"max equation residual={trace['maximum_residual']:.3e}; "
          f"min rank-three row volume={trace['minimum_row_volume']:.3e}")
    print("crossings: T2=100, pillowcase=50; "
          f"iota error<={crossing_diagnostics['worst_match']:.3e}; "
          f"min angle sine={crossing_diagnostics['minimum_angle_sine']:.3e}")
    print(f"deck census={dict(sorted(EXPECTED_DECK_CENSUS.items()))}")
    print("analytic deck-(2,2) nodes:")
    for point in deck_two:
        print(f"  P=({point[0]:.7f},{point[1]:.7f})")
    print("KWZ: 74 T2 hits -> 37 dots -> 3 cancellations -> 31/30; "
          "delta^2=0")
    print("root subdivisions 20/40/80 agree; reversed interval relabeling "
          f"intertwines all {len(base_pre['delta'])} arrows")
    print("affine iota quotient: worst height gap="
          f"{max(item['maximum_pair_gap'] for item in quotient_diagnostics):.3e}; "
          "worst relative gap="
          f"{max(item['maximum_relative_gap'] for item in quotient_diagnostics):.3e}")
    return {
        "deck_two": deck_two,
        "signature": signatures[0],
        "matching": matching,
        "curve": curve,
        "lift": lift,
        "switch_data": switch_data,
        "switch_matching": switch_matching,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--perturbation", type=float, action="append",
                        help="positive C3 perturbation; may be repeated")
    parser.add_argument("--step", type=float, default=0.04,
                        help="maximum pseudo-arclength step (default: 0.04)")
    parser.add_argument("--stability", action="store_true",
                        help="also run the independent t=0.02 trace")
    parser.add_argument(
        "--switch-census", action="store_true",
        help="smooth and pair both exact closures at all 50 analytic nodes")
    args = parser.parse_args(argv)
    perturbations = args.perturbation or [0.015]
    if args.stability and 0.02 not in perturbations:
        perturbations.append(0.02)
    if any(value <= 0.0 for value in perturbations):
        raise ValueError("perturbations must be positive")

    base = base_precurve()
    names = pl_deck_two_nodes(base[0], base[1])
    print("PL connector--main nodes with short deck class (2,2): "
          + ", ".join(names))
    results = [run_one(value, args.step, base) for value in perturbations]
    switch_censuses = []
    if args.switch_census:
        switch_censuses = [
            analytic_switch_census(result, base) for result in results
        ]
    if len(results) > 1:
        if not all(result["signature"] == results[0]["signature"]
                   for result in results[1:]):
            raise AssertionError("analytic type-D object changed with perturbation")
        print("\n[PASS] all requested perturbations give the same type-D object")
        if args.switch_census:
            if not all(census["signature"] == switch_censuses[0]["signature"]
                       for census in switch_censuses[1:]):
                raise AssertionError(
                    "analytic switch census changed with perturbation")
            print("[PASS] all requested perturbations give the same "
                  "singleton-switch and two-closure census")

    print("\nESTABLISHED BY THIS REGRESSION:")
    print("  the corrected floating-point C3 trace and the q=7 PL curve encode "
          "strictly isomorphic 31-generator type-D base objects")
    print("  the analytic trace has exactly two deck-(2,2) nodes, corresponding "
          "to the uniquely characterized PL nodes S18/S25")
    if args.switch_census:
        print("  among all 43 module-preserving singleton smoothings of the "
              "actual curve, the rank pair (9,31) selects only S18/S25")
        print("  those switches are independent, and all 2^43 F2 sums "
              "satisfy the type-D Maurer--Cartan equation")
    print("NOT ESTABLISHED:")
    print("  interval-rigorous existence/completeness of the numerical trace")
    print("  the localized corrections to Smith's typeset C3 formulas as an "
          "author-verified erratum")
    print("  which switch sum (if any) is selected by rigid C3 figure-eight "
          "counts, or that presentation-changing outputs are excluded")
    print("  the CHKK instanton--pillowcase assignment or its bounding-cochain "
          "support")
    return 0


if __name__ == "__main__":
    sys.exit(main())
