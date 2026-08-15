#!/usr/bin/env python3
"""Exact singular-limit root and action certificate for Paper III.

Paper III is the unconditional action paper in ``atiyah_floer/paper3``.  The
separate full Gate-A project is preserved in ``atiyah_floer/research_program``.

The calculation is deliberately independent of the floating continuation in
``q7_action_audit.py``.  It uses ``fractions.Fraction`` throughout to

* reconstruct the first-arc ``A union H_0 union H_pi`` path word;
* isolate the algebraic roots of the seam-normal displacement equations;
* reduce 58 exact singular overlaps to the 50 transverse roots that persist;
* determine the degree-one ordering and short deck class; and
* evaluate the 25 half-period representatives by the area--residue formula.

All action values printed by this module are coefficients of ``pi^2``.  The
only transcendental enclosures used to locate the normal roots are rational
bounds for pi and cosine.  They are constructed by Machin's formula and
Taylor bounds, rather than by converting floating-point output to fractions.

Run from the repository root with either

    python3 pillowcase/q7_exact_actions.py

or

    ~/git/data_platform/.venv/bin/python pillowcase/q7_action_audit.py \
        --exact-limit
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial


Interval = tuple[F, F]
Point = tuple[F, F]


# ---------------------------------------------------------------------------
# Rational interval arithmetic
# ---------------------------------------------------------------------------


def i_point(value: int | F) -> Interval:
    value = F(value)
    return value, value


def i_add(first: Interval, second: Interval) -> Interval:
    return first[0] + second[0], first[1] + second[1]


def i_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def i_sub(first: Interval, second: Interval) -> Interval:
    return i_add(first, i_neg(second))


def i_mul(first: Interval, second: Interval) -> Interval:
    products = (
        first[0] * second[0], first[0] * second[1],
        first[1] * second[0], first[1] * second[1],
    )
    return min(products), max(products)


def i_inv(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        raise ZeroDivisionError(f"interval contains zero: {value}")
    reciprocals = (1 / value[0], 1 / value[1])
    return min(reciprocals), max(reciprocals)


def i_div(first: Interval, second: Interval) -> Interval:
    return i_mul(first, i_inv(second))


def i_scale(coefficient: int | F, value: Interval) -> Interval:
    return i_mul(i_point(coefficient), value)


def i_square(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        return F(0), max(value[0] * value[0], value[1] * value[1])
    return i_mul(value, value)


def strict_sign(value: Interval) -> int:
    if value[1] < 0:
        return -1
    if value[0] > 0:
        return 1
    raise AssertionError(f"interval does not certify a strict sign: {value}")


def arctan_reciprocal_bounds(denominator: int, terms: int) -> Interval:
    """Alternating rational bounds for arctan(1/denominator)."""
    total = F(0)
    for index in range(terms):
        term = F(
            1,
            (2 * index + 1) * denominator ** (2 * index + 1),
        )
        total += term if index % 2 == 0 else -term
    omitted = F(
        1,
        (2 * terms + 1) * denominator ** (2 * terms + 1),
    )
    if terms % 2 == 0:
        return total, total + omitted
    return total - omitted, total


def pi_bounds() -> Interval:
    """Machin: pi/4 = 4 arctan(1/5) - arctan(1/239)."""
    first = arctan_reciprocal_bounds(5, 24)
    second = arctan_reciprocal_bounds(239, 8)
    return i_scale(4, i_sub(i_scale(4, first), second))


PI = pi_bounds()


def cosine_point_bounds(argument: F, terms: int = 36) -> Interval:
    """Rational Taylor enclosure for cos(argument), with a tail bound."""
    total = F(0)
    square = argument * argument
    power = F(1)
    for index in range(terms):
        if index:
            power *= square
        term = power / factorial(2 * index)
        total += term if index % 2 == 0 else -term
    omitted = abs(argument) ** (2 * terms) / factorial(2 * terms)
    return total - omitted, total + omitted


def cosine_pi_multiple_bounds(multiplier: F) -> Interval:
    """Enclose cos(multiplier*pi) for 0 <= multiplier <= 1."""
    if not 0 <= multiplier <= 1:
        raise ValueError(f"cosine monotonicity range violated: {multiplier}")
    angle = i_scale(multiplier, PI)
    at_upper = cosine_point_bounds(angle[1])
    at_lower = cosine_point_bounds(angle[0])
    return at_upper[0], at_lower[1]


# The unique root 2 cos(pi/7) of u^3-u^2-2u+1 in this interval.
U: Interval = (F(1801937735, 10**9), F(1801937737, 10**9))


def minimal_polynomial(value: F) -> F:
    return value**3 - value**2 - 2 * value + 1


def polynomial_coefficients(name: str) -> tuple[Interval, Interval, Interval]:
    """Coefficients A,B,C of the three squared normal equations."""
    u_squared = i_square(U)
    if name == "12":
        return (
            i_sub(i_sub(u_squared, i_scale(2, U)), i_point(1)),
            i_add(
                i_add(i_scale(F(-5, 2), u_squared), i_scale(F(3, 2), U)),
                i_point(F(5, 2)),
            ),
            i_sub(
                i_sub(i_scale(F(1, 4), u_squared), i_scale(F(1, 2), U)),
                i_point(F(1, 4)),
            ),
        )
    if name == "13":
        return (
            i_sub(i_sub(u_squared, i_scale(F(1, 2), U)), i_point(F(3, 2))),
            i_add(
                i_add(i_scale(F(-3, 4), u_squared), i_scale(F(-1, 2), U)),
                i_point(2),
            ),
            i_sub(
                i_sub(i_scale(F(1, 4), u_squared), i_scale(F(1, 8), U)),
                i_point(F(3, 8)),
            ),
        )
    if name == "23":
        return (
            i_sub(i_add(u_squared, U), i_point(2)),
            i_sub(i_sub(u_squared, i_scale(F(5, 2), U)), i_point(2)),
            i_sub(
                i_add(i_scale(F(1, 4), u_squared), i_scale(F(1, 4), U)),
                i_point(F(1, 2)),
            ),
        )
    raise KeyError(name)


def polynomial_value(name: str, value: Interval) -> Interval:
    quadratic, linear, constant = polynomial_coefficients(name)
    return i_add(
        i_add(i_mul(quadratic, i_square(value)), i_mul(linear, value)),
        constant,
    )


def polynomial_derivative(name: str, value: Interval) -> Interval:
    quadratic, linear, _ = polynomial_coefficients(name)
    return i_add(i_scale(2, i_mul(quadratic, value)), linear)


@dataclass(frozen=True)
class RootBox:
    polynomial: str | None
    cosine: Interval
    scaled_angle: Interval


ROOT_BOXES: dict[str, RootBox] = {
    "A": RootBox(
        "13",
        (F(178447, 10**6), F(178449, 10**6)),
        (F(9300, 1000), F(9301, 1000)),
    ),
    "B": RootBox(
        "12",
        (F(-123491, 10**6), F(-123489, 10**6)),
        (F(11327, 1000), F(11328, 1000)),
    ),
    "C": RootBox(
        None,
        (F(-311746, 10**6), F(-311743, 10**6)),
        (F(12619, 1000), F(12620, 1000)),
    ),
    "D": RootBox(
        "23",
        (F(722519, 10**6), F(722522, 10**6)),
        (F(5102, 1000), F(5103, 1000)),
    ),
    "E": RootBox(
        "23",
        (F(346009, 10**6), F(346012, 10**6)),
        (F(8138, 1000), F(8139, 1000)),
    ),
    "F": RootBox(
        None,
        (F(111259, 10**6), F(111262, 10**6)),
        (F(9754, 1000), F(9755, 1000)),
    ),
    "G": RootBox(
        None,
        (F(450483, 10**6), F(450486, 10**6)),
        (F(7376, 1000), F(7377, 1000)),
    ),
}


def check_angle_box(cosine: Interval, scaled_angle: Interval) -> None:
    """Prove arccos(cosine)*21/pi lies in scaled_angle."""
    cosine_at_lower = cosine_pi_multiple_bounds(scaled_angle[0] / 21)
    cosine_at_upper = cosine_pi_multiple_bounds(scaled_angle[1] / 21)
    if not (
        cosine_at_upper[1] < cosine[0]
        and cosine[1] < cosine_at_lower[0]
    ):
        raise AssertionError(
            (cosine, scaled_angle, cosine_at_lower, cosine_at_upper)
        )


def row_centers() -> dict[int, Interval]:
    """Return a_i=cos(theta_1)cos(theta_2) for the six seam rows."""
    u_squared = i_square(U)
    first = i_scale(F(-1, 4), i_sub(u_squared, i_point(2)))
    second = i_scale(
        F(1, 4), i_sub(i_sub(u_squared, U), i_point(1))
    )
    third = i_scale(F(1, 4), U)
    return {1: first, 2: second, 3: third,
            4: third, 5: second, 6: first}


ROW_CENTER = row_centers()
ROW_TYPE = {1: 1, 6: 1, 2: 2, 5: 2, 3: 3, 4: 3}
# The traversed lift has cos(alpha)>0 on rows 1--3 and <0 on rows 4--6.
ROW_ALPHA_SIGN = {1: 1, 2: 1, 3: 1, 4: -1, 5: -1, 6: -1}
ZERO_ROOT_BY_TYPE = {1: "C", 2: "F", 3: "G"}
SQUARED_ROOTS_BY_TYPES = {
    (1, 2): ("B",),
    (1, 3): ("A",),
    (2, 3): ("D", "E"),
}


def row_width_squared(row: int) -> Interval:
    """b_i^2=sin^2(theta_1)sin^2(theta_2)."""
    center = ROW_CENTER[row]
    return i_scale(F(3, 4), i_sub(i_point(1), i_scale(4, i_square(center))))


def normal_sign(row: int, cosine: Interval) -> int:
    """Sign of G_i=2(a_i-cos(theta))/cos(alpha_i)."""
    numerator = i_sub(ROW_CENTER[row], cosine)
    return ROW_ALPHA_SIGN[row] * strict_sign(numerator)


def normal_roots(first_row: int, second_row: int) -> tuple[str, ...]:
    """Physical simple roots of G_first-G_second on 0<theta<pi."""
    first_type = ROW_TYPE[first_row]
    second_type = ROW_TYPE[second_row]
    if first_type == second_type:
        if ROW_ALPHA_SIGN[first_row] == ROW_ALPHA_SIGN[second_row]:
            raise AssertionError("coincident normal functions are not isolated")
        return (ZERO_ROOT_BY_TYPE[first_type],)
    pair = tuple(sorted((first_type, second_type)))
    candidates = SQUARED_ROOTS_BY_TYPES[pair]
    return tuple(
        label for label in candidates
        if normal_sign(first_row, ROOT_BOXES[label].cosine)
        == normal_sign(second_row, ROOT_BOXES[label].cosine)
    )


def check_algebraic_root_certificate() -> None:
    """Run every rational enclosure used by the seam-root census."""
    if not minimal_polynomial(U[0]) < 0 < minimal_polynomial(U[1]):
        raise AssertionError("u interval does not bracket the algebraic root")
    derivative = i_sub(
        i_sub(i_scale(3, i_square(U)), i_scale(2, U)), i_point(2)
    )
    if derivative[0] <= 0:
        raise AssertionError("minimal polynomial is not monotone on u interval")
    trig_u = i_scale(2, cosine_pi_multiple_bounds(F(1, 7)))
    if not U[0] < trig_u[0] <= trig_u[1] < U[1]:
        raise AssertionError("2 cos(pi/7) is not enclosed by U")

    for label in ("A", "B", "D", "E"):
        box = ROOT_BOXES[label]
        left = polynomial_value(box.polynomial, i_point(box.cosine[0]))
        right = polynomial_value(box.polynomial, i_point(box.cosine[1]))
        if not ((left[1] < 0 < right[0]) or (right[1] < 0 < left[0])):
            raise AssertionError((label, left, right))
        strict_sign(polynomial_derivative(box.polynomial, box.cosine))
        check_angle_box(box.cosine, box.scaled_angle)

    # P12 and P13 have one physical root; Vieta places the other outside [-1,1].
    for polynomial, label, side in (("12", "B", -1), ("13", "A", 1)):
        quadratic, linear, _ = polynomial_coefficients(polynomial)
        root_sum = i_div(i_neg(linear), quadratic)
        other = i_sub(root_sum, ROOT_BOXES[label].cosine)
        if side < 0 and other[1] >= -1:
            raise AssertionError((polynomial, other))
        if side > 0 and other[0] <= 1:
            raise AssertionError((polynomial, other))

    # P23 is quadratic and the two disjoint sign-changing boxes exhaust it.
    if ROOT_BOXES["E"].cosine[1] >= ROOT_BOXES["D"].cosine[0]:
        raise AssertionError("P23 root boxes are not disjoint")

    # The three opposite-sign row pairs meet exactly where their common
    # numerator vanishes: x=a_1,a_2,a_3.
    for row, label in ((1, "C"), (2, "F"), (3, "G")):
        center = ROW_CENTER[row]
        box = ROOT_BOXES[label]
        if not box.cosine[0] < center[0] <= center[1] < box.cosine[1]:
            raise AssertionError((row, label, center, box.cosine))
        check_angle_box(box.cosine, box.scaled_angle)
        if row_width_squared(row)[0] <= 0:
            raise AssertionError("zero normal root lies at a fold")


# ---------------------------------------------------------------------------
# Exact singular first arc and its overlaps
# ---------------------------------------------------------------------------


# Heights are n=21*theta/pi.  A missing key means the diagonal A sheet simply
# crosses the seam.  Values are the signed vertical H_0/H_pi jumps.
JUMPS = {
    1: -6, 4: -12, 5: 6, 7: -14, 9: 12, 10: 14,
    11: 14, 12: 12, 14: -14, 16: 6, 17: -12, 20: -6,
}
ROW_BY_SEAM = {
    1: 4, 4: 1, 5: 4, 7: 5, 9: 6, 10: 2,
    11: 5, 12: 1, 14: 2, 16: 3, 17: 6, 20: 3,
}


def singular_path() -> tuple[list[Point], dict[int, tuple[F, F]]]:
    """Build the first arc in coordinates (gamma/pi, theta/pi)."""
    vertices: list[Point] = [(F(0), F(0))]
    intervals: dict[int, tuple[F, F]] = {0: (F(0), F(0))}
    height = F(0)
    for seam in range(1, 22):
        before = height + 10
        vertices.append((F(seam), before / 21))
        after = before + JUMPS.get(seam, 0)
        intervals[seam] = tuple(sorted((before, after)))
        if after != before:
            vertices.append((F(seam), after / 21))
        height = after
    if vertices[-1] != (F(21), F(10)) or sum(JUMPS.values()) != 0:
        raise AssertionError("singular path does not have deck class (21,10)")
    return vertices, intervals


VERTICES, SEAM_INTERVALS = singular_path()


@dataclass(frozen=True)
class Overlap:
    first: int
    second: int
    offset: int
    lower: F
    upper: F
    root: str | None
    root_height: Interval


def transformed_angle_boxes(label: str, lower: F, upper: F) -> list[Interval]:
    """All translates/reflections of a root-angle box contained in a range."""
    base = ROOT_BOXES[label].scaled_angle
    found: list[Interval] = []
    for period in range(-6, 8):
        candidates = (
            i_add(base, i_point(42 * period)),
            i_add(i_neg(base), i_point(42 * (period + 1))),
        )
        for candidate in candidates:
            if lower < candidate[0] and candidate[1] < upper:
                found.append(candidate)
    return found


def all_singular_overlaps() -> list[tuple[int, int, int, F, F]]:
    overlaps = []
    for first in range(1, 21):
        for second in range(first + 1, 21):
            for offset in range(-210, 211, 42):
                lower = max(
                    SEAM_INTERVALS[first][0],
                    SEAM_INTERVALS[second][0] - offset,
                )
                upper = min(
                    SEAM_INTERVALS[first][1],
                    SEAM_INTERVALS[second][1] - offset,
                )
                if lower <= upper:
                    overlaps.append((first, second, offset, lower, upper))
    return overlaps


def persistent_overlaps() -> tuple[list[Overlap], list[tuple[int, int]]]:
    """Classify the 58 singular overlaps by exact local transversality."""
    persistent: list[Overlap] = []
    excluded: list[tuple[int, int]] = []
    for first, second, offset, lower, upper in all_singular_overlaps():
        first_vertical = first in JUMPS
        second_vertical = second in JUMPS
        if first_vertical != second_vertical:
            # The A sheet meets the interior of one vertical circle.
            if lower != upper:
                raise AssertionError("point--interval overlap is not a point")
            vertical = first if first_vertical else second
            vertical_height = lower if first_vertical else lower + offset
            interval = SEAM_INTERVALS[vertical]
            if min(vertical_height - interval[0],
                   interval[1] - vertical_height) < 2:
                raise AssertionError("point--interval root approaches a fold")
            persistent.append(Overlap(
                first, second, offset, lower, upper, None,
                (lower, lower),
            ))
            continue
        if not (first_vertical and second_vertical):
            raise AssertionError("unexpected A--A singular overlap")

        labels = normal_roots(ROW_BY_SEAM[first], ROW_BY_SEAM[second])
        located = [
            (label, boxes[0])
            for label in labels
            if len(boxes := transformed_angle_boxes(label, lower, upper)) == 1
        ]
        if not located:
            excluded.append((first, second))
            continue
        if len(located) != 1:
            raise AssertionError((first, second, located))
        label, height = located[0]
        persistent.append(Overlap(
            first, second, offset, lower, upper, label, height
        ))

    if len(all_singular_overlaps()) != 58:
        raise AssertionError("singular overlap census changed")
    if len(persistent) != 50 or len(excluded) != 8:
        raise AssertionError((len(persistent), len(excluded)))
    expected_excluded = {
        (1, 17), (4, 14), (4, 20), (5, 9),
        (7, 17), (9, 11), (10, 12), (12, 16),
    }
    if set(excluded) != expected_excluded:
        raise AssertionError((excluded, expected_excluded))
    return persistent, excluded


# ---------------------------------------------------------------------------
# Degree-one orderings and deck classes
# ---------------------------------------------------------------------------


def reduced_sine_sign(scaled_angle: Interval) -> int:
    quotient = scaled_angle[0].numerator // scaled_angle[0].denominator
    period = quotient // 42
    reduced = i_sub(scaled_angle, i_point(42 * period))
    if 0 < reduced[0] and reduced[1] < 21:
        return 1
    if 21 < reduced[0] and reduced[1] < 42:
        return -1
    raise AssertionError(f"angle box meets a zero of sine: {scaled_angle}")


def logarithmic_normal_term(row: int, cosine: Interval) -> Interval:
    """The row-dependent part of d/dtheta log|G_i|."""
    difference = i_sub(ROW_CENTER[row], cosine)
    denominator = i_sub(
        row_width_squared(row),
        i_square(i_sub(cosine, ROW_CENTER[row])),
    )
    if denominator[0] <= 0:
        raise AssertionError("normal root reaches a fold")
    return i_add(i_inv(difference), i_div(difference, denominator))


def degree_one_order(overlap: Overlap) -> str:
    """Determine the Akaho--Joyce degree-one branch ordering exactly."""
    first_vertical = overlap.first in JUMPS
    second_vertical = overlap.second in JUMPS
    if first_vertical != second_vertical:
        # The diagonal A direction is (1,10/21).  The circle direction has
        # signed vertical component JUMPS[seam].
        determinant = (
            -JUMPS[overlap.first]
            if first_vertical else JUMPS[overlap.second]
        )
        return "A->B" if determinant > 0 else "B->A"

    first_row = ROW_BY_SEAM[overlap.first]
    second_row = ROW_BY_SEAM[overlap.second]
    sine_sign = reduced_sine_sign(overlap.root_height)
    if overlap.root in ("C", "F", "G"):
        # Here G_second=-G_first and both vanish simply.
        derivative_sign = sine_sign * ROW_ALPHA_SIGN[first_row]
    else:
        cosine = ROOT_BOXES[overlap.root].cosine
        g_sign = normal_sign(first_row, cosine)
        log_difference = i_sub(
            logarithmic_normal_term(first_row, cosine),
            logarithmic_normal_term(second_row, cosine),
        )
        derivative_sign = g_sign * sine_sign * strict_sign(log_difference)
    determinant_sign = (
        (1 if JUMPS[overlap.first] > 0 else -1)
        * (1 if JUMPS[overlap.second] > 0 else -1)
        * derivative_sign
    )
    return "A->B" if determinant_sign > 0 else "B->A"


def canonical_deck(vector: tuple[int, int]) -> tuple[int, int]:
    if vector[0] < 0 or (vector[0] == 0 and vector[1] < 0):
        return -vector[0], -vector[1]
    return vector


def short_deck(overlap: Overlap) -> tuple[int, int]:
    cross_arc = (overlap.second - overlap.first) % 2 == 1
    if cross_arc:
        forward = (
            (overlap.second + 21 - overlap.first) // 2,
            overlap.offset // 42 + 5,
        )
    else:
        forward = (
            (overlap.second - overlap.first) // 2,
            overlap.offset // 42,
        )
    reverse = (21 - forward[0], 10 - forward[1])
    chosen = (
        forward
        if sum(map(abs, forward)) <= sum(map(abs, reverse))
        else reverse
    )
    return canonical_deck(chosen)


# ---------------------------------------------------------------------------
# Exact area--residue evaluation
# ---------------------------------------------------------------------------


def path_position(point: Point) -> F:
    seam, height = point
    indices = [
        index for index, vertex in enumerate(VERTICES)
        if vertex[0] == seam
    ]
    if len(indices) == 1:
        if height != VERTICES[indices[0]][1]:
            raise AssertionError((point, VERTICES[indices[0]]))
        return F(indices[0])
    first, second = indices
    start = VERTICES[first][1]
    end = VERTICES[second][1]
    fraction = (height - start) / (end - start)
    if not 0 <= fraction <= 1:
        raise AssertionError((point, start, end))
    return F(first) + fraction


def forward_subpath(start: Point, end: Point) -> list[Point]:
    start_position = path_position(start)
    end_position = path_position(end)
    if start_position > end_position:
        return list(reversed(forward_subpath(end, start)))
    path = [start]
    first_vertex = start_position.numerator // start_position.denominator + 1
    last_vertex = end_position.numerator // end_position.denominator
    for index in range(first_vertex, last_vertex + 1):
        vertex = VERTICES[index]
        if vertex != path[-1] and vertex != end:
            path.append(vertex)
    if end != path[-1]:
        path.append(end)
    return path


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def point_on_segment(point: Point, start: Point, end: Point) -> bool:
    cross = (
        (end[0] - start[0]) * (point[1] - start[1])
        - (end[1] - start[1]) * (point[0] - start[0])
    )
    return (
        cross == 0
        and min(start[0], end[0]) <= point[0] <= max(start[0], end[0])
        and min(start[1], end[1]) <= point[1] <= max(start[1], end[1])
    )


def assert_no_lattice_points_on_polygon(polygon: list[Point]) -> None:
    """Check genericity by scanning each edge's own integer bounding box."""
    for start, end in zip(polygon, polygon[1:]):
        for x in range(floor_fraction(min(start[0], end[0])),
                       ceil_fraction(max(start[0], end[0])) + 1):
            for y in range(floor_fraction(min(start[1], end[1])),
                           ceil_fraction(max(start[1], end[1])) + 1):
                point = F(x), F(y)
                if point_on_segment(point, start, end):
                    raise AssertionError(
                        "auxiliary polygon meets a deleted lattice point: "
                        f"{point}"
                    )


def winding_number(polygon: list[Point], point: Point) -> int:
    px, py = point
    winding = 0
    for first, second in zip(polygon, polygon[1:]):
        x1, y1 = first
        x2, y2 = second
        cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1)
        if y1 <= py < y2 and cross > 0:
            winding += 1
        elif y2 <= py < y1 and cross < 0:
            winding -= 1
    return winding


def normalized_loop_integral(
    path: list[Point],
    closure_offset: Point = (F(137, 1000), F(271, 1000)),
) -> tuple[F, F, int]:
    """Return (integral/pi^2, signed area, total lattice winding)."""
    start = path[0]
    end = path[-1]
    translate = end[0] - start[0], end[1] - start[1]
    if any(value.denominator != 1 or value % 2 for value in translate):
        raise AssertionError(f"path does not close on the 2pi torus: {translate}")
    reflected_start = -start[0] - translate[0], -start[1] - translate[1]
    midpoint = (
        (reflected_start[0] + start[0]) / 2 + closure_offset[0],
        (reflected_start[1] + start[1]) / 2 + closure_offset[1],
    )
    polygon = list(path)
    polygon.extend([
        (midpoint[0] + translate[0], midpoint[1] + translate[1]),
        (-start[0], -start[1]),
    ])
    polygon.extend([(-point[0], -point[1]) for point in path[1:]])
    polygon.extend([midpoint, start])
    area = F(1, 2) * sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(polygon, polygon[1:])
    )

    minimum_x = min(point[0] for point in polygon)
    maximum_x = max(point[0] for point in polygon)
    minimum_y = min(point[1] for point in polygon)
    maximum_y = max(point[1] for point in polygon)
    lattice = [
        (F(x), F(y))
        for x in range(floor_fraction(minimum_x) - 1,
                       ceil_fraction(maximum_x) + 2)
        for y in range(floor_fraction(minimum_y) - 1,
                       ceil_fraction(maximum_y) + 2)
    ]
    assert_no_lattice_points_on_polygon(polygon)
    winding = sum(winding_number(polygon, point) for point in lattice)
    return F(1, 2) * (area - winding), area, winding


def exact_integral(start: Point, end: Point) -> tuple[F, F, int, F]:
    """Integral from start to end, plus its area/residue certificate.

    The last return value is the multiplier applied to the symmetrized-loop
    formula: one for an ordinary torus loop and one half for a half-shift loop.
    """
    path = forward_subpath(start, end)
    delta = end[0] - start[0], end[1] - start[1]
    if delta[0].denominator != 1 or delta[1].denominator != 1:
        raise AssertionError(f"nonintegral endpoint shift: {delta}")
    if delta[1] % 2:
        raise AssertionError(f"theta shift is not a torus period: {delta}")
    if delta[0] % 2 == 0:
        multiplier = F(1)
        working_path = path
    else:
        multiplier = F(1, 2)
        working_path = list(path)
        working_path.extend([
            (point[0] + delta[0], point[1] + delta[1])
            for point in path[1:]
        ])
    primary = normalized_loop_integral(working_path)
    check = normalized_loop_integral(
        working_path, (F(319, 1000), F(-183, 1000))
    )
    if primary[0] != check[0]:
        raise AssertionError((primary, check))
    return multiplier * primary[0], primary[1], primary[2], multiplier


@dataclass(frozen=True)
class ActionRow:
    overlap: Overlap
    mate: tuple[int, int]
    deck: tuple[int, int]
    ordering: str
    area: F
    winding: int
    loop_multiplier: F
    action: F


def half_period_mate(overlap: Overlap) -> tuple[int, int]:
    return 21 - overlap.second, 21 - overlap.first


def action_for_overlap(overlap: Overlap) -> tuple[F, F, int, F]:
    height = (overlap.lower + overlap.upper) / 2
    start = F(overlap.first), height / 21
    end = F(overlap.second), (height + overlap.offset) / 21
    integral, area, winding, multiplier = exact_integral(start, end)
    ordering = degree_one_order(overlap)
    action = -integral if ordering == "A->B" else integral
    return action, area, winding, multiplier


EXPECTED_CENSUS = Counter({
    F(-13, 14): 2, F(-6, 7): 2, F(-11, 14): 2,
    F(-5, 7): 2, F(-9, 14): 2, F(-13, 21): 2,
    F(-1, 2): 6, F(-10, 21): 2, F(-19, 42): 2,
    F(-3, 7): 2, F(-5, 14): 2, F(-1, 3): 2,
    F(-13, 42): 2, F(-3, 14): 2, F(-4, 21): 2,
    F(-1, 6): 2, F(-1, 14): 2, F(-1, 21): 2,
    F(-1, 42): 2, F(5, 42): 2, F(3, 14): 2,
    F(5, 14): 2, F(3, 7): 2,
})


def exact_action_rows() -> tuple[list[ActionRow], Counter]:
    persistent, _ = persistent_overlaps()
    by_pair = {(item.first, item.second): item for item in persistent}
    visited: set[tuple[int, int]] = set()
    rows: list[ActionRow] = []
    census: Counter = Counter()
    for overlap in persistent:
        pair = overlap.first, overlap.second
        if pair in visited:
            continue
        mate_pair = half_period_mate(overlap)
        mate = by_pair.get(mate_pair)
        if mate is None or mate.offset != overlap.offset:
            raise AssertionError((pair, mate_pair, overlap.offset))
        visited.update((pair, mate_pair))
        action, area, winding, multiplier = action_for_overlap(overlap)
        mate_action, _, _, _ = action_for_overlap(mate)
        if action != mate_action:
            raise AssertionError((pair, mate_pair, action, mate_action))

        rows.append(ActionRow(
            overlap=overlap,
            mate=mate_pair,
            deck=short_deck(overlap),
            ordering=degree_one_order(overlap),
            area=area,
            winding=winding,
            loop_multiplier=multiplier,
            action=action,
        ))
        census[action] += 2
    if len(rows) != 25 or len(visited) != 50:
        raise AssertionError((len(rows), len(visited)))
    if census != EXPECTED_CENSUS:
        raise AssertionError((census, EXPECTED_CENSUS))
    return rows, census


def fraction_text(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def print_report() -> None:
    check_algebraic_root_certificate()
    persistent, excluded = persistent_overlaps()
    rows, census = exact_action_rows()
    point_interval = sum(
        (item.first in JUMPS) != (item.second in JUMPS)
        for item in persistent
    )
    interval_interval = len(persistent) - point_interval

    print("EXACT q=7 SINGULAR-LIMIT ROOT/ACTION CERTIFICATE")
    print(f"  rational pi enclosure width: {float(PI[1] - PI[0]):.3e}")
    print(f"  u=2*cos(pi/7) enclosure: [{float(U[0]):.10f}, "
          f"{float(U[1]):.10f}]")
    print("  normal-root boxes (n=21*theta/pi):")
    for label, box in ROOT_BOXES.items():
        print(
            f"    {label}: x in [{float(box.cosine[0]):.6f},"
            f"{float(box.cosine[1]):.6f}], n in "
            f"[{float(box.scaled_angle[0]):.3f},"
            f"{float(box.scaled_angle[1]):.3f}]"
        )
    print("  first-arc vertical jumps (seam:signed height): "
          + ", ".join(f"{seam}:{jump:+d}" for seam, jump in JUMPS.items()))
    print(
        f"  singular overlaps=58; persistent={len(persistent)} "
        f"({point_interval} point--interval + "
        f"{interval_interval} interval--interval); excluded={len(excluded)}"
    )
    print("  excluded seam pairs: "
          + ", ".join(f"({first},{second})" for first, second in excluded))
    print("\n  25 half-period representatives; each row has multiplicity two")
    print("  pair   mate   off root deck   deg   area      wind fac action/pi^2")
    for row in rows:
        item = row.overlap
        root = row.overlap.root or "pt"
        print(
            f"  ({item.first:2d},{item.second:2d}) "
            f"({row.mate[0]:2d},{row.mate[1]:2d}) "
            f"{item.offset:+4d}  {root:>2s}  "
            f"({row.deck[0]:2d},{row.deck[1]:1d}) "
            f"{row.ordering:>4s} "
            f"{fraction_text(row.area):>9s} "
            f"{row.winding:5d} "
            f"{fraction_text(row.loop_multiplier):>3s} "
            f"{fraction_text(row.action):>10s}"
        )
    negative = sum(count for value, count in census.items() if value < 0)
    positive = sum(count for value, count in census.items() if value > 0)
    zero = sum(count for value, count in census.items() if value == 0)
    print("\n  exact action census:")
    for value, count in sorted(census.items()):
        print(f"    {fraction_text(value):>6s} * pi^2 : {count}")
    print(f"  signs: negative={negative}, positive={positive}, zero={zero}")
    print("  minimum absolute limiting action: (1/42)*pi^2")
    print("PASS: exact 42-negative/8-positive singular-limit census")


def main() -> int:
    print_report()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
