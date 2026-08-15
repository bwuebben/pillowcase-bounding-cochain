#!/usr/bin/env python3
"""Smith's holonomy-perturbed C3 correspondence in quaternion coordinates.

This is an equation-level model, not the regular-homotopy connector model in
``resolve.py``.  It implements the Wirtinger presentation and boundary labels
in Smith, arXiv:2412.06066v1, Proposition 4.1 and Figure 20.

There are source-level typographical inconsistencies in the preprint.
Proposition 4.1 and Equation (4.1.1) use ``lambda = a c^{-1} x``, but the
expansion immediately below has the opposite sign.  In Figure 20, the S2
top-right meridian is printed ``x p^{-1} x c x p x^{-1}``; replacing the
second ``x`` after ``c`` by ``x^{-1}`` makes it the conjugate
``(x p^{-1} x)c(x p^{-1} x)^{-1}``, makes it reduce to ``c`` at ``t=0``, and
makes the redundant Wirtinger relation hold exactly.  ``phi_smith_printed``
retains Equation (4.1.5) as a diagnostic; the geometric model here is defined
by the group presentation with this localized boundary-word correction.

For parameters ``(gamma, theta, alpha, beta)`` set

    a = i,  b = exp(gamma k)i,  c = exp(theta k)i,
    x = sin(alpha)cos(beta)i + sin(alpha)sin(beta)j + cos(alpha)k.

The perturbation and remaining meridians are

    p = exp(t Im(a c^{-1} x)),
    d = c a^{-1} b,
    y = x p^{-1} a^{-1} p b.

The perturbed traceless variety is the scalar equation Re(y)=0.  Figure 20
gives the ordered boundary quadruples (a_i,b_i,c_i,d_i):

    S1 = (p^{-1}ap, b, x, y),
    S2 = (xp^{-1}xpx^{-1}, y, xp^{-1}xcx^{-1}px^{-1}, d),
    S3 = (a, b, c, d).

The corrected S2 word also equals the meridian reconstructed from the other
three by the exact pillowcase relation.  The displayed word instead equals
its negative in every traceless SU(2) representation and reduces to ``-c`` at
``t=0``.

The pillowcase convention is Smith's: gamma_i is the smaller angle from a_i
to b_i, and theta_i is the angle from a_i to c_i oriented toward b_i.  This
agrees with the repository convention that Q_{p/q} is
``q*theta - p*gamma = 0 (mod 2*pi)``.

Run the equation/regression checks with

    python3 pillowcase/c3_perturbed.py
"""
from __future__ import annotations

import math
import sys


PI = math.pi
TAU = 2.0 * PI
I = (0.0, 1.0, 0.0, 0.0)


def qmul(first, second):
    """Hamilton product of quaternions in (real,i,j,k) coordinates."""
    a, b, c, d = first
    e, f, g, h = second
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


def qprod(*factors):
    result = (1.0, 0.0, 0.0, 0.0)
    for factor in factors:
        result = qmul(result, factor)
    return result


def qconj(quaternion):
    return (quaternion[0], -quaternion[1], -quaternion[2], -quaternion[3])


def qdistance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def _cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _norm(vector):
    return math.sqrt(_dot(vector, vector))


def _unit(vector):
    norm = _norm(vector)
    if norm == 0.0:
        raise ValueError("zero vector has no direction")
    return tuple(value / norm for value in vector)


def _clamp(value):
    return max(-1.0, min(1.0, value))


def plane_meridian(angle):
    """The pure unit quaternion exp(angle*k)i."""
    return (0.0, math.cos(angle), math.sin(angle), 0.0)


def sphere_meridian(alpha, beta):
    return (
        0.0,
        math.sin(alpha) * math.cos(beta),
        math.sin(alpha) * math.sin(beta),
        math.cos(alpha),
    )


def exp_imaginary(quaternion, scale):
    """exp(scale * Im(quaternion)), including the zero-vector limit."""
    vector = quaternion[1:]
    norm = _norm(vector)
    if norm < 1.0e-15:
        return (1.0, 0.0, 0.0, 0.0)
    angle = scale * norm
    factor = math.sin(angle) / norm
    return (math.cos(angle), *(factor * value for value in vector))


def representation(gamma, theta, alpha, beta, perturbation):
    """Return Smith's seven named meridian holonomies."""
    a = I
    b = plane_meridian(gamma)
    c = plane_meridian(theta)
    x = sphere_meridian(alpha, beta)
    d = qprod(c, qconj(a), b)
    # Proposition 4.1 and Equation (4.1.1): lambda = a c^{-1} x.
    lam = qprod(a, qconj(c), x)
    p = exp_imaginary(lam, perturbation)
    y = qprod(x, qconj(p), qconj(a), p, b)
    return {"a": a, "b": b, "c": c, "d": d,
            "x": x, "y": y, "p": p, "lambda": lam}


def phi_direct(gamma, theta, alpha, beta, perturbation):
    """Phi_t = Re(y); its zero set is the perturbed traceless variety."""
    return representation(gamma, theta, alpha, beta, perturbation)["y"][0]


def _sinc(value):
    return 1.0 if abs(value) < 1.0e-12 else math.sin(value) / value


def phi_smith_printed(gamma, theta, alpha, beta, perturbation):
    """The preprint's printed Equation (4.1.5), retained as a diagnostic."""
    sa, ca = math.sin(alpha), math.cos(alpha)
    st, ct = math.sin(theta), math.cos(theta)
    sb, cb = math.sin(beta), math.cos(beta)
    n = perturbation * math.sqrt(max(0.0, 1.0 - ca * ca * st * st))
    cn = math.cos(n)
    sn = _sinc(n)
    f_value = (
        2.0 * cn * sn * perturbation
        * (sa * sa * sb * math.sin(theta - beta) - ca * ca * ct)
        + 2.0 * sn * sn * perturbation * perturbation
        * sa * sa * cb * ca * st * math.cos(beta - theta)
    )
    g_value = (
        ca * cn * cn
        - 2.0 * cn * sn * perturbation
        * sa * sa * cb * math.sin(theta - beta)
        + sn * sn * perturbation * perturbation
        * ca * sa * sa * math.cos(beta - theta) * math.cos(beta + theta)
    )
    return math.cos(gamma) * f_value + math.sin(gamma) * g_value


def boundary_quadruples(rep):
    """The repaired ordered (a_i,b_i,c_i,d_i) boundary tuples.

    Figure 20 gives S2's first meridian ``x p^{-1} x p x^{-1}``.  Its third
    word is missing an inverse on the ``x`` immediately following ``c``.
    With that correction it is a conjugate of ``c`` and agrees identically
    with the meridian reconstructed from the boundary relation.
    """
    a, b, c, d = (rep[name] for name in ("a", "b", "c", "d"))
    x, y, p = (rep[name] for name in ("x", "y", "p"))
    s1_a = qprod(qconj(p), a, p)
    s2_a = qprod(x, qconj(p), x, p, qconj(x))
    s2_c = qprod(x, qconj(p), x, c, qconj(x), p, qconj(x))
    return {
        "S1": (s1_a, b, x, y),
        "S2": (s2_a, y, s2_c, d),
        "S3": (a, b, c, d),
    }


def printed_s2_third_meridian(rep):
    """The inconsistent top-right S2 word displayed in Figure 20."""
    c, x, p = (rep[name] for name in ("c", "x", "p"))
    return qprod(x, qconj(p), x, c, x, p, qconj(x))


def pillowcase_relation(quadruple):
    """Residual norm for Smith's pillowcase relation a_i c_i^{-1}=b_i d_i^{-1}."""
    a, b, c, d = quadruple
    return qdistance(qprod(a, qconj(c)), qprod(b, qconj(d)))


def pillowcase_coordinates(quadruple, tolerance=2.0e-8):
    """Read Smith's (gamma,theta) coordinates from a traceless boundary tuple.

    The returned representative has gamma in [0,pi] and theta in [0,2pi).
    At gamma=0 or pi, theta is folded to [0,pi], as required by the pillowcase
    involution.
    """
    if pillowcase_relation(quadruple) > 20.0 * tolerance:
        raise ValueError("boundary tuple violates the pillowcase relation")
    if any(abs(quaternion[0]) > tolerance for quaternion in quadruple):
        raise ValueError("boundary tuple is not traceless")
    a, b, c, _ = (quaternion[1:] for quaternion in quadruple)
    gamma = math.acos(_clamp(_dot(a, b)))
    normal = _cross(a, b)
    if _norm(normal) <= tolerance:
        theta = math.acos(_clamp(_dot(a, c)))
    else:
        normal = _unit(normal)
        theta = math.atan2(_dot(normal, _cross(a, c)), _dot(a, c)) % TAU
    return gamma, theta


def pillowcase_coordinates_relaxed(quadruple, tolerance=1.0e-12):
    """A smooth local extension of pillowcase coordinates off Phi_t=0.

    Newton iteration evaluates nearby tuples for which y is not quite
    traceless.  Normalizing the imaginary parts gives a valid extension near
    a solution away from gamma=0,pi.  No mathematical claim is made about the
    off-zero-set values themselves.
    """
    axes = []
    for quaternion in quadruple[:3]:
        vector = quaternion[1:]
        if _norm(vector) <= tolerance:
            raise ValueError("nearby boundary meridian has no imaginary direction")
        axes.append(_unit(vector))
    a, b, c = axes
    gamma = math.acos(_clamp(_dot(a, b)))
    normal = _cross(a, b)
    if _norm(normal) <= tolerance:
        raise ValueError("relaxed pillowcase chart reached a gamma seam")
    normal = _unit(normal)
    theta = math.atan2(_dot(normal, _cross(a, c)), _dot(a, c)) % TAU
    return gamma, theta


def _unwrap_near(value, target):
    return value + TAU * round((target - value) / TAU)


def lifted_coordinates_near(quadruple, target):
    """Choose the T^2 lift of a boundary point nearest ``target``.

    The two lifts are (gamma,theta) and (-gamma,-theta).  Integer translates
    are included when comparing to the requested unwrapped branch.
    """
    gamma, theta = pillowcase_coordinates_relaxed(quadruple)
    candidates = ((gamma, theta), ((-gamma) % TAU, (-theta) % TAU))
    lifted = [
        (_unwrap_near(first, target[0]), _unwrap_near(second, target[1]))
        for first, second in candidates
    ]
    return min(
        lifted,
        key=lambda point: ((point[0] - target[0]) ** 2
                           + (point[1] - target[1]) ** 2),
    )


def boundary_coordinates(gamma, theta, alpha, beta, perturbation):
    rep = representation(gamma, theta, alpha, beta, perturbation)
    if abs(rep["y"][0]) > 2.0e-8:
        raise ValueError("parameters do not lie on Phi_t=0")
    return {
        name: pillowcase_coordinates(quadruple)
        for name, quadruple in boundary_quadruples(rep).items()
    }


def rational_residual(coordinates, numerator, denominator):
    """Smooth residual for q*theta-p*gamma = 0 modulo 2*pi."""
    gamma, theta = coordinates
    return math.sin(0.5 * (denominator * theta - numerator * gamma))


def composition_residuals(parameters, perturbation, first_denominator=3,
                          second_denominator=7, sheets=(0, 0), targets=None):
    """Three lifted equations cutting out Q_{1/q1} C3-composed with Q_{1/q2}.

    ``parameters`` are (gamma,theta,alpha,beta).  The sheet integers select
    ``q_i*theta_i-gamma_i=2*pi*k_i``.  ``targets`` select continuous T^2
    lifts during Newton iteration.
    """
    gamma, theta, alpha, beta = parameters
    phi = phi_direct(gamma, theta, alpha, beta, perturbation)
    if targets is None:
        targets = (
            (gamma, (gamma + TAU * sheets[0]) / first_denominator),
            (gamma, (gamma + TAU * sheets[1]) / second_denominator),
        )
    rep = representation(gamma, theta, alpha, beta, perturbation)
    boundaries = boundary_quadruples(rep)
    first = lifted_coordinates_near(boundaries["S1"], targets[0])
    second = lifted_coordinates_near(boundaries["S2"], targets[1])
    return (
        phi,
        first_denominator * first[1] - first[0] - TAU * sheets[0],
        second_denominator * second[1] - second[0] - TAU * sheets[1],
    )


def _solve_linear3(matrix, vector, pivot_tolerance=1.0e-12):
    """Solve a real 3x3 system by partial-pivot Gaussian elimination."""
    augmented = [list(row) + [value] for row, value in zip(matrix, vector)]
    minimum_pivot = math.inf
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(augmented[row][column]))
        size = abs(augmented[pivot][column])
        minimum_pivot = min(minimum_pivot, size)
        if size <= pivot_tolerance:
            raise ValueError("singular 3x3 Newton matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(3):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return tuple(augmented[row][3] for row in range(3)), minimum_pivot


def _solve_linear(matrix, vector, pivot_tolerance=1.0e-12):
    """Solve a small dense real system by partial-pivot elimination."""
    size = len(vector)
    if len(matrix) != size or any(len(row) != size for row in matrix):
        raise ValueError("linear system is not square")
    augmented = [list(row) + [value] for row, value in zip(matrix, vector)]
    minimum_pivot = math.inf
    for column in range(size):
        pivot = max(range(column, size),
                    key=lambda row: abs(augmented[row][column]))
        magnitude = abs(augmented[pivot][column])
        minimum_pivot = min(minimum_pivot, magnitude)
        if magnitude <= pivot_tolerance:
            raise ValueError("singular Newton matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return tuple(augmented[row][size] for row in range(size)), minimum_pivot


def _determinant3(matrix):
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _jacobian(function, point, output_size):
    matrix = [[0.0] * len(point) for _ in range(output_size)]
    for column in range(len(point)):
        step = 1.5e-6 * max(1.0, abs(point[column]))
        upper = list(point)
        lower = list(point)
        upper[column] += step
        lower[column] -= step
        r_upper = function(upper)
        r_lower = function(lower)
        for row in range(output_size):
            matrix[row][column] = (r_upper[row] - r_lower[row]) / (2.0 * step)
    return matrix


def _null_tangent3x4(jacobian):
    """The oriented cofactor vector spanning the nullspace of a rank-3 3x4 matrix."""
    tangent = []
    for removed in range(4):
        minor = [[row[column] for column in range(4) if column != removed]
                 for row in jacobian]
        tangent.append(((-1.0) ** removed) * _determinant3(minor))
    norm = math.sqrt(sum(value * value for value in tangent))
    if norm <= 1.0e-12:
        raise ValueError("composition Jacobian lost rank")
    return tuple(value / norm for value in tangent)


def unwrapped_composition_residuals(parameters, perturbation,
                                    denominators=(3, 7)):
    """Global lifted equations along the full q1*q2-sheet main component."""
    gamma, theta, alpha, beta = parameters
    rep = representation(gamma, theta, alpha, beta, perturbation)
    boundaries = boundary_quadruples(rep)
    first_target = (gamma, gamma / denominators[0])
    second_target = (gamma, gamma / denominators[1])
    first = lifted_coordinates_near(boundaries["S1"], first_target)
    second = lifted_coordinates_near(boundaries["S2"], second_target)
    return (
        rep["y"][0],
        denominators[0] * first[1] - first[0],
        denominators[1] * second[1] - second[0],
    )


def unperturbed_unwrapped(gamma, denominators=(3, 7)):
    first = gamma / denominators[0]
    second = gamma / denominators[1]
    return (gamma, first + second, PI / 2.0, first)


def qpower(quaternion, exponent):
    if exponent < 0:
        return qpower(qconj(quaternion), -exponent)
    result = (1.0, 0.0, 0.0, 0.0)
    factor = quaternion
    power = exponent
    while power:
        if power & 1:
            result = qmul(result, factor)
        factor = qmul(factor, factor)
        power >>= 1
    return result


def rational_quaternion_residual(quadruple, denominator):
    """Quaternion residual for Q_{1/q}: (a^{-1}c)^q=a^{-1}b."""
    a, b, c, _ = quadruple
    left = qpower(qprod(qconj(a), c), denominator)
    right = qprod(qconj(a), b)
    return tuple(first - second for first, second in zip(left, right))


def algebraic_composition_residuals(parameters, perturbation,
                                    denominators=(3, 7)):
    """Chart-free corrected-C3 plus two rational-tangle equations."""
    rep = representation(*parameters, perturbation)
    boundaries = boundary_quadruples(rep)
    return (
        rep["y"][0],
        *rational_quaternion_residual(boundaries["S1"], denominators[0]),
        *rational_quaternion_residual(boundaries["S2"], denominators[1]),
    )


def solve_algebraic_at_gamma(gamma, perturbation, initial=None,
                             denominators=(3, 7), tolerance=2.0e-10):
    """Overdetermined fixed-gamma initializer using chart-free equations."""
    base = unperturbed_unwrapped(gamma, denominators)
    values = list(initial[1:] if initial is not None else base[1:])

    def residual(reduced):
        return algebraic_composition_residuals(
            (gamma, *reduced), perturbation, denominators)

    for _ in range(30):
        current = residual(values)
        norm = math.sqrt(sum(value * value for value in current))
        if norm <= tolerance:
            return (gamma, *values)
        jacobian = _jacobian(residual, values, len(current))
        normal = [[sum(jacobian[row][i] * jacobian[row][j]
                       for row in range(len(current)))
                   for j in range(3)] for i in range(3)]
        right = tuple(-sum(jacobian[row][i] * current[row]
                           for row in range(len(current)))
                      for i in range(3))
        delta, _ = _solve_linear(normal, right, pivot_tolerance=1.0e-14)
        scale = 1.0
        for _ in range(14):
            candidate = [value + scale * change
                         for value, change in zip(values, delta)]
            candidate_norm = math.sqrt(sum(
                value * value for value in residual(candidate)))
            if candidate_norm < norm:
                values = candidate
                break
            scale *= 0.5
        else:
            raise ValueError(
                f"algebraic fixed-gamma initializer stalled at {norm:.3e}")
    raise ValueError("algebraic fixed-gamma initializer did not converge")


def _maximal_row_tangent(jacobian):
    """Null tangent from the best-conditioned independent three-row minor."""
    best_tangent = None
    best_norm = 0.0
    row_count = len(jacobian)
    for first in range(row_count - 2):
        for second in range(first + 1, row_count - 1):
            for third in range(second + 1, row_count):
                rows = [jacobian[first], jacobian[second], jacobian[third]]
                cofactors = []
                for removed in range(4):
                    minor = [[row[column] for column in range(4)
                              if column != removed] for row in rows]
                    cofactors.append(((-1.0) ** removed) * _determinant3(minor))
                norm = math.sqrt(sum(value * value for value in cofactors))
                if norm > best_norm:
                    best_norm = norm
                    best_tangent = tuple(value / norm for value in cofactors)
    if best_tangent is None or best_norm <= 1.0e-10:
        raise ValueError("algebraic composition Jacobian lost rank three")
    return best_tangent, best_norm


def _arclength_correct_overdetermined(function, prediction, tangent,
                                      tolerance=2.0e-10):
    point = list(prediction)
    minimum_pivot = math.inf
    for iteration in range(15):
        geometric = function(point)
        hyperplane = sum((point[index] - prediction[index]) * tangent[index]
                         for index in range(4))
        norm = math.sqrt(sum(value * value for value in geometric)
                         + hyperplane * hyperplane)
        if norm <= tolerance:
            return tuple(point), iteration, minimum_pivot, norm
        jacobian = _jacobian(function, point, len(geometric))
        normal = [[sum(jacobian[row][i] * jacobian[row][j]
                       for row in range(len(geometric)))
                   + tangent[i] * tangent[j]
                   for j in range(4)] for i in range(4)]
        right = tuple(
            -sum(jacobian[row][i] * geometric[row]
                 for row in range(len(geometric)))
            - tangent[i] * hyperplane
            for i in range(4)
        )
        delta, pivot = _solve_linear(normal, right, pivot_tolerance=1.0e-15)
        minimum_pivot = min(minimum_pivot, pivot)
        scale = 1.0
        for _ in range(12):
            candidate = [value + scale * change
                         for value, change in zip(point, delta)]
            candidate_geometric = function(candidate)
            candidate_hyperplane = sum(
                (candidate[index] - prediction[index]) * tangent[index]
                for index in range(4))
            candidate_norm = math.sqrt(
                sum(value * value for value in candidate_geometric)
                + candidate_hyperplane * candidate_hyperplane)
            if candidate_norm < norm:
                point = candidate
                break
            scale *= 0.5
        else:
            raise ValueError("overdetermined pseudo-arclength corrector stalled")
    raise ValueError("overdetermined pseudo-arclength corrector did not converge")


def solve_unwrapped_at_gamma(gamma, perturbation, initial=None,
                             denominators=(3, 7), tolerance=2.0e-11):
    """Fixed-gamma initialization for pseudo-arclength continuation."""
    base = unperturbed_unwrapped(gamma, denominators)
    values = list(initial[1:] if initial is not None else base[1:])

    def residual(reduced):
        return unwrapped_composition_residuals(
            (gamma, *reduced), perturbation, denominators)

    for iteration in range(25):
        current = residual(values)
        norm = math.sqrt(sum(value * value for value in current))
        if norm <= tolerance:
            return (gamma, *values)
        jacobian = _jacobian(residual, values, 3)
        delta, _ = _solve_linear(jacobian, tuple(-value for value in current))
        scale = 1.0
        for _ in range(14):
            candidate = [value + scale * change
                         for value, change in zip(values, delta)]
            candidate_norm = math.sqrt(sum(
                value * value for value in residual(candidate)))
            if candidate_norm < norm:
                values = candidate
                break
            scale *= 0.5
        else:
            raise ValueError("fixed-gamma initializer stalled")
    raise ValueError("fixed-gamma initializer did not converge")


def _arclength_correct(function, prediction, tangent, tolerance=2.0e-11):
    point = list(prediction)
    minimum_pivot = math.inf
    for iteration in range(12):
        geometric = function(point)
        hyperplane = sum((point[index] - prediction[index]) * tangent[index]
                         for index in range(4))
        residual = tuple(geometric) + (hyperplane,)
        norm = math.sqrt(sum(value * value for value in residual))
        if norm <= tolerance:
            return tuple(point), iteration, minimum_pivot, norm
        jacobian = _jacobian(function, point, 3) + [list(tangent)]
        delta, pivot = _solve_linear(jacobian, tuple(-value for value in residual))
        minimum_pivot = min(minimum_pivot, pivot)
        scale = 1.0
        for _ in range(10):
            candidate = [value + scale * change
                         for value, change in zip(point, delta)]
            candidate_geometric = function(candidate)
            candidate_hyperplane = sum(
                (candidate[index] - prediction[index]) * tangent[index]
                for index in range(4))
            candidate_norm = math.sqrt(
                sum(value * value for value in candidate_geometric)
                + candidate_hyperplane * candidate_hyperplane)
            if candidate_norm < norm:
                point = candidate
                break
            scale *= 0.5
        else:
            raise ValueError("pseudo-arclength corrector stalled")
    raise ValueError("pseudo-arclength corrector did not converge")


def trace_corrected_main_arclength(perturbation, step=0.035,
                                   denominators=(3, 7), max_steps=12000):
    """Chart-free pseudo-arclength trace of the repaired main component."""
    period = math.lcm(*denominators) * TAU
    start_gamma = 0.13
    start = solve_algebraic_at_gamma(start_gamma, perturbation,
                                     denominators=denominators)
    function = lambda point: algebraic_composition_residuals(
        point, perturbation, denominators)
    tangent, row_volume = _maximal_row_tangent(
        _jacobian(function, start, len(function(start))))
    if tangent[0] < 0.0:
        tangent = tuple(-value for value in tangent)

    points = [start]
    diagnostics = {
        "maximum_residual": 0.0,
        "minimum_corrector_pivot": math.inf,
        "maximum_corrector_iterations": 0,
        "minimum_gamma_tangent": tangent[0],
        "maximum_gamma_tangent": tangent[0],
        "minimum_row_volume": row_volume,
        "maximum_tangent_residual": 0.0,
    }
    target_gamma = start_gamma + period
    current_step = step
    for _ in range(max_steps):
        prediction = tuple(points[-1][index] + current_step * tangent[index]
                           for index in range(4))
        try:
            corrected, iterations, pivot, residual = (
                _arclength_correct_overdetermined(
                    function, prediction, tangent))
        except ValueError:
            current_step *= 0.5
            if current_step < step / 32.0:
                raise ValueError(
                    "pseudo-arclength corrector failed near "
                    f"(gamma,theta,alpha,beta)={points[-1]}"
                )
            continue
        points.append(corrected)
        jacobian = _jacobian(function, corrected, len(function(corrected)))
        new_tangent, row_volume = _maximal_row_tangent(jacobian)
        if sum(a * b for a, b in zip(new_tangent, tangent)) < 0.0:
            new_tangent = tuple(-value for value in new_tangent)
        tangent = new_tangent
        diagnostics["maximum_residual"] = max(
            diagnostics["maximum_residual"], residual)
        diagnostics["minimum_corrector_pivot"] = min(
            diagnostics["minimum_corrector_pivot"], pivot)
        diagnostics["maximum_corrector_iterations"] = max(
            diagnostics["maximum_corrector_iterations"], iterations)
        diagnostics["minimum_gamma_tangent"] = min(
            diagnostics["minimum_gamma_tangent"], tangent[0])
        diagnostics["maximum_gamma_tangent"] = max(
            diagnostics["maximum_gamma_tangent"], tangent[0])
        diagnostics["minimum_row_volume"] = min(
            diagnostics["minimum_row_volume"], row_volume)
        diagnostics["maximum_tangent_residual"] = max(
            diagnostics["maximum_tangent_residual"],
            max(abs(sum(row[column] * tangent[column]
                        for column in range(4))) for row in jacobian),
        )
        if iterations <= 3:
            current_step = min(step, current_step * 1.25)
        if corrected[0] >= target_gamma:
            break
    else:
        raise ValueError(
            "pseudo-arclength trace exceeded its step budget near "
            f"(gamma,theta,alpha,beta)={points[-1]} with tangent={tangent}"
        )

    expected_translate = (
        period,
        period * (1.0 / denominators[0] + 1.0 / denominators[1]),
        0.0,
        period / denominators[0],
    )
    closure_error = math.sqrt(sum(
        (points[-1][index] - points[0][index] - expected_translate[index]) ** 2
        for index in range(4)))
    diagnostics["raw_closure_error"] = closure_error
    if closure_error > max(0.05, 4.0 * step):
        raise ValueError(
            "pseudo-arclength trace did not return to the expected deck "
            f"translate (error={closure_error:.3e})"
        )
    points[-1] = tuple(points[0][index] + expected_translate[index]
                       for index in range(4))
    diagnostics["closure_overshoot_error"] = 0.0
    diagnostics["steps"] = len(points) - 1
    torus_curve = [(point[0] % TAU, point[1] % TAU) for point in points]
    torus_curve.append(torus_curve[0])
    return torus_curve, points, diagnostics


def unperturbed_branch(gamma, sheets, denominators=(3, 7)):
    """Exact t=0 point on one lifted rational sheet."""
    first = (gamma + TAU * sheets[0]) / denominators[0]
    second = (gamma + TAU * sheets[1]) / denominators[1]
    return (gamma, first + second, PI / 2.0, first)


def solve_branch_at_gamma(gamma, perturbation, sheets, initial=None,
                          denominators=(3, 7), tolerance=2.0e-11,
                          max_iterations=20):
    """Newton solve one corrected-C3 branch with output gamma held fixed."""
    base = unperturbed_branch(gamma, sheets, denominators)
    parameters = list(initial[1:] if initial is not None else base[1:])
    targets = (
        (gamma, base[3]),
        (gamma, base[1] - base[3]),
    )

    def residual(values):
        return composition_residuals(
            (gamma, *values), perturbation,
            first_denominator=denominators[0],
            second_denominator=denominators[1],
            sheets=sheets, targets=targets,
        )

    minimum_pivot = math.inf
    for iteration in range(max_iterations):
        current = residual(parameters)
        norm = math.sqrt(sum(value * value for value in current))
        if norm <= tolerance:
            return {
                "parameters": (gamma, *parameters),
                "residual": norm,
                "iterations": iteration,
                "minimum_pivot": minimum_pivot,
            }
        jacobian = [[0.0] * 3 for _ in range(3)]
        for column in range(3):
            step = 2.0e-6 * max(1.0, abs(parameters[column]))
            upper = list(parameters)
            lower = list(parameters)
            upper[column] += step
            lower[column] -= step
            r_upper = residual(upper)
            r_lower = residual(lower)
            for row in range(3):
                jacobian[row][column] = (r_upper[row] - r_lower[row]) / (2.0 * step)
        delta, pivot = _solve_linear3(jacobian, tuple(-value for value in current))
        minimum_pivot = min(minimum_pivot, pivot)
        delta_norm = math.sqrt(sum(value * value for value in delta))
        if delta_norm > 0.75:
            delta = tuple(value * 0.75 / delta_norm for value in delta)

        accepted = False
        scale = 1.0
        for _ in range(12):
            candidate = [value + scale * change
                         for value, change in zip(parameters, delta)]
            candidate_residual = residual(candidate)
            candidate_norm = math.sqrt(
                sum(value * value for value in candidate_residual))
            if candidate_norm < norm:
                parameters = candidate
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            raise ValueError(
                f"Newton line search stalled at residual {norm:.3e}")
    raise ValueError(f"Newton did not converge after {max_iterations} iterations")


def trace_corrected_main(perturbation, samples_per_sheet=120,
                         denominators=(3, 7)):
    """Trace the corrected generic main component on the pillowcase T^2 cover.

    Holding the output gamma fixed cuts the component into q1*q2 sheets.  At
    gamma=2*pi the sheet (k1,k2) joins (k1+1,k2+1) at gamma=0.  For coprime
    denominators this permutation is one cycle.  A half-grid offset avoids the
    coordinate singularities at gamma=0 and pi.
    """
    first_denominator, second_denominator = denominators
    expected_length = math.lcm(first_denominator, second_denominator)
    if expected_length != first_denominator * second_denominator:
        raise ValueError("trace_corrected_main currently expects coprime denominators")
    gammas = [TAU * (index + 0.5) / samples_per_sheet
              for index in range(samples_per_sheet)]
    points = []
    diagnostics = {
        "maximum_residual": 0.0,
        "minimum_pivot": math.inf,
        "maximum_iterations": 0,
    }
    sheets = (0, 0)
    seen = []
    for _ in range(expected_length):
        if sheets in seen:
            raise AssertionError("sheet endpoint permutation closed too early")
        seen.append(sheets)
        previous = None
        for gamma in gammas:
            try:
                solution = solve_branch_at_gamma(
                    gamma, perturbation, sheets, initial=previous,
                    denominators=denominators,
                )
            except ValueError:
                # A fixed lifted chart can jump across a seam.  Restart from
                # the exact t=0 point before declaring that the branch is lost.
                solution = solve_branch_at_gamma(
                    gamma, perturbation, sheets, initial=None,
                    denominators=denominators,
                )
            previous = solution["parameters"]
            _, theta, _, _ = solution["parameters"]
            points.append((gamma % TAU, theta % TAU))
            diagnostics["maximum_residual"] = max(
                diagnostics["maximum_residual"], solution["residual"])
            diagnostics["minimum_pivot"] = min(
                diagnostics["minimum_pivot"], solution["minimum_pivot"])
            diagnostics["maximum_iterations"] = max(
                diagnostics["maximum_iterations"], solution["iterations"])
        sheets = ((sheets[0] + 1) % first_denominator,
                  (sheets[1] + 1) % second_denominator)
    if sheets != (0, 0) or len(seen) != expected_length:
        raise AssertionError("sheet endpoint permutation did not form the expected cycle")
    points.append(points[0])
    diagnostics["sheet_cycle"] = seen
    diagnostics["points"] = len(points) - 1
    return points, diagnostics


def partial_correspondence_triple(perturbation, tolerance=1.0e-14):
    """Return the three exact-branch parameters in the partial triple fiber.

    Fixing ``gamma=0`` and ``theta=pi/2``, the choices ``beta=0`` and
    ``beta=pi`` make the corrected C3 traceless equation automatic.  On the
    torus cover the Q_{1/3} equation becomes, respectively,

        3*(pi/2-alpha) + 4*t*sin(alpha) = 0,
        3*(pi/2+alpha) - 4*t*sin(alpha) = 2*pi*k,

    with k=1,2 in the second line.  Their t=0 roots pi/2, pi/6,
    and 5*pi/6 are simple, so Newton continuation is unambiguous for the
    small perturbations used here.
    """
    if not 0.0 <= perturbation < 0.25:
        raise ValueError("triple-point continuation expects 0 <= t < 1/4")
    specifications = (
        (0.0, PI / 2.0, 0),
        (PI, PI / 6.0, 1),
        (PI, 5.0 * PI / 6.0, 2),
    )
    points = []
    for beta, center, sheet in specifications:
        alpha = center
        for _ in range(20):
            if beta == 0.0:
                value = (3.0 * (PI / 2.0 - alpha)
                         + 4.0 * perturbation * math.sin(alpha))
                derivative = -3.0 + 4.0 * perturbation * math.cos(alpha)
            else:
                value = (3.0 * (PI / 2.0 + alpha)
                         - 4.0 * perturbation * math.sin(alpha)
                         - TAU * sheet)
                derivative = 3.0 - 4.0 * perturbation * math.cos(alpha)
            alpha -= value / derivative
            if abs(value) <= tolerance:
                break
        else:
            raise ValueError("partial triple-point scalar solve did not converge")
        points.append((0.0, PI / 2.0, alpha, beta))
    return points


def _angle_distance(first, second):
    difference = (first - second) % TAU
    return min(difference, TAU - difference)


def run_checks():
    failures = []

    def check(label, error, tolerance=1.0e-9):
        passed = error <= tolerance
        print(f"[{'PASS' if passed else 'FAIL'}] {label}: error={error:.3e}")
        if not passed:
            failures.append((label, error, tolerance))

    samples = [
        (0.37, 1.11, 0.42, 2.03, 0.07),
        (1.28, 5.17, 1.39, 0.83, 0.13),
        (2.71, 2.22, 2.41, 4.91, 0.19),
        (5.43, 3.74, 0.91, 5.72, 0.23),
    ]
    printed_gap = max(
        abs(phi_direct(*sample) - phi_smith_printed(*sample))
        for sample in samples
    )
    print("[SOURCE DISCREPANCY] exact Wirtinger Re(y) versus printed "
          f"Equation (4.1.5): max gap={printed_gap:.3e}")

    longitude_error = 0.0
    commutator_error = 0.0
    for gamma, theta, alpha, beta, perturbation in samples:
        rep = representation(gamma, theta, alpha, beta, perturbation)
        expected = (
            math.cos(alpha) * math.sin(theta),
            math.sin(alpha) * math.cos(theta - beta),
            -math.sin(alpha) * math.sin(theta - beta),
            math.cos(alpha) * math.cos(theta),
        )
        longitude_error = max(longitude_error, qdistance(rep["lambda"], expected))
        commutator_error = max(
            commutator_error,
            qdistance(qmul(rep["p"], rep["lambda"]),
                      qmul(rep["lambda"], rep["p"])),
        )
    check("lambda=a*c^{-1}*x has the independently expanded coordinates",
          longitude_error, 3.0e-15)
    check("the perturbation meridian p commutes with its longitude",
          commutator_error, 3.0e-15)

    relation_errors = {"S1": 0.0, "S2": 0.0, "S3": 0.0}
    corrected_word_error = 0.0
    displayed_word_sign_error = 0.0
    for sample in samples:
        rep = representation(*sample)
        boundaries = boundary_quadruples(rep)
        for name, quadruple in boundaries.items():
            relation_errors[name] = max(
                relation_errors[name], pillowcase_relation(quadruple))
        s2_a, _, s2_c, _ = boundaries["S2"]
        reconstructed = qprod(rep["d"], qconj(rep["y"]), s2_a)
        corrected_word_error = max(
            corrected_word_error, qdistance(s2_c, reconstructed))
        displayed = printed_s2_third_meridian(rep)
        displayed_word_sign_error = max(
            displayed_word_sign_error,
            qdistance(displayed, tuple(-value for value in s2_c)),
        )
    for name in ("S1", "S2", "S3"):
        check(f"Figure-20 {name} tuple satisfies its boundary relation",
              relation_errors[name], 3.0e-15)
    check("corrected S2 word equals the relation-derived meridian",
          corrected_word_error, 3.0e-15)
    check("displayed S2 word is the negative of the corrected meridian",
          displayed_word_sign_error, 3.0e-15)

    coordinate_error = 0.0
    for gamma, theta, beta in (
            (0.31, 1.77, 0.42),
            (1.19, 5.21, 2.31),
            (2.73, 3.07, 5.62)):
        coordinates = boundary_coordinates(gamma, theta, PI / 2.0, beta, 0.0)
        expected = {
            "S1": (gamma, beta % TAU),
            "S2": (gamma, (theta - beta) % TAU),
            "S3": (gamma, theta % TAU),
        }
        for name in coordinates:
            coordinate_error = max(
                coordinate_error,
                abs(coordinates[name][0] - expected[name][0]),
                _angle_distance(coordinates[name][1], expected[name][1]),
            )
    check("t=0 boundary coordinates recover theta-addition", coordinate_error)

    rational_error = 0.0
    for gamma in (0.23, 0.91, 1.83, 2.87):
        for first_sheet in range(3):
            for second_sheet in range(7):
                beta = (gamma + TAU * first_sheet) / 3.0
                second_theta = (gamma + TAU * second_sheet) / 7.0
                theta = beta + second_theta
                coordinates = boundary_coordinates(
                    gamma, theta, PI / 2.0, beta, 0.0)
                rational_error = max(
                    rational_error,
                    abs(rational_residual(coordinates["S1"], 1, 3)),
                    abs(rational_residual(coordinates["S2"], 1, 7)),
                    abs(math.sin(0.5 * (21.0 * coordinates["S3"][1]
                                        - 10.0 * coordinates["S3"][0]))),
                )
    check("Q_1/3 and Q_1/7 compose to the slope-10/21 subtorus",
          rational_error, 2.0e-14)

    algebraic_sheet_error = 0.0
    for gamma in (0.23, 2.17, 5.81, 8.29):
        parameters = unperturbed_unwrapped(gamma)
        algebraic_sheet_error = max(
            algebraic_sheet_error,
            max(abs(value) for value in
                algebraic_composition_residuals(parameters, 0.0)),
        )
    check("chart-free quaternion equations recover the t=0 main lift",
          algebraic_sheet_error, 3.0e-14)

    initialized = solve_algebraic_at_gamma(0.37, 0.015)
    initialized_error = max(abs(value) for value in
                            algebraic_composition_residuals(initialized, 0.015))
    check("nonzero-t chart-free algebraic initializer", initialized_error,
          3.0e-10)

    triple_parameters = partial_correspondence_triple(0.015)
    triple_c3_error = 0.0
    triple_rational_error = 0.0
    triple_target_error = 0.0
    triple_p1_gammas = []
    triple_x = []
    for parameters in triple_parameters:
        rep = representation(*parameters, 0.015)
        boundaries = boundary_quadruples(rep)
        coordinates = {
            name: pillowcase_coordinates(quadruple)
            for name, quadruple in boundaries.items()
        }
        triple_c3_error = max(triple_c3_error, abs(rep["y"][0]))
        triple_rational_error = max(
            triple_rational_error,
            max(abs(value) for value in
                rational_quaternion_residual(boundaries["S1"], 3)),
        )
        for name in ("S2", "S3"):
            triple_target_error = max(
                triple_target_error,
                abs(coordinates[name][0]),
                _angle_distance(coordinates[name][1], PI / 2.0),
            )
        triple_p1_gammas.append(coordinates["S1"][0])
        triple_x.append(rep["x"])
    triple_domain_separation = min(
        qdistance(triple_x[first], triple_x[second])
        for first in range(3) for second in range(first + 1, 3)
    )
    check("partial triple fiber lies on corrected C3", triple_c3_error,
          3.0e-14)
    check("partial triple fiber lies on Q_1/3", triple_rational_error,
          3.0e-13)
    check("three partial-correspondence branches have target "
          "((0,pi/2),(0,pi/2))", triple_target_error, 3.0e-13)
    check("partial triple preimages are distinct and avoid P1 corners",
          0.0 if (triple_domain_separation > 0.5
                  and min(triple_p1_gammas) > 1.0e-3) else 1.0)

    print(f"\nC3 equation checks: {len(failures)} failure(s)")
    if failures:
        for label, error, tolerance in failures:
            print(f"  {label}: {error} > {tolerance}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(run_checks())
