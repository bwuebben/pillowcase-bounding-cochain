#!/usr/bin/env python3
"""
tangles.py -- the compass-pinned tangle calculus on the pillowcase (RESEARCH_LOG sec 23).

Resolves the sec-22 blocker (puncture labeling) by DERIVATION, not convention-guessing:

  * A traceless SU(2) rep of the 4-punctured sphere forces the four meridian axes to be
    COPLANAR (Re(X1 X2 X3) = 0 is the scalar triple product of the axes), so every rep has
    the binary-dihedral normal form beta = (0, gamma, theta+gamma, theta) of sec 1.
  * The product relation x1 x2 x3 x4 = 1 gives exact angle identities:
    West-pair angle = East-pair angle, North = South.
  * SLOT CONVENTION: (x1,x2,x3,x4) = (NW, SW, SE, NE), counterclockwise from NW.
    gamma = angle(x1,x2) = WEST angle; theta = angle(x1,x4) = NORTH angle.
  * Conway sum T1+T2 glues T1's East punctures to T2's West punctures. East=West inside
    each tangle => the glued angle is SHARED: the sum is the fiber product over gamma,
    and rotating T2's frame onto the glued pair shows theta ADDS:
        L(T1+T2) = {(g, t1+t2) : (g,t1) in L1~, (g,t2) in L2~}  on the torus double cover,
    where L~ is the FULL iota-invariant lift (both lifts are legitimate gluing frames).
  * Twist calculus (derived by the reflection recursion b -> 2*phi - b):
    East twists (slots SE,NE)  : line theta = n*gamma   (Conway fraction n, slope n)
    South twists (slots SW,SE) : line theta = gamma/n   (fraction 1/n, slope 1/n)
    rotate = cyclic slot shift : (gamma,theta) -> (-theta,gamma), i.e. F -> -1/F.
    Fraction p/q <-> slope p/q <-> class (q,p); numerator closure pairs with theta=0
    (det = p), denominator closure with gamma=0 (det = q).  pcase.py Layer 1's
    "(p,q) paired with the infinity line" is the 90-degree-rotated presentation of
    the same dictionary -- both det-correct, hence the old steep/shallow ambiguity.

Sec 24 adds the SEAM FIBER-CIRCLE rule (the sum components the generic fiber product
misses): when both curves touch a seam gamma in {0,pi} at non-corner points, the S^1
stabilizer of the aligned glued axes gives a circle component mapping 2:1 onto the
seam arc [|th1-th2|, pi-|pi-(th1+th2)|], with
    cos th_out(alpha) = cos th1 cos th2 + sin th1 sin th2 cos alpha
-- derived here, validated by an explicit-SO(3) 3D engine (glue3d_*), and identical
to Smith's published formula (arXiv:2412.06066).

Everything here is grounded rep theory (quaternion words), no Smith conventions.
Validation battery in __main__; run `python3 tangles.py` (50 checks).
"""
import math
from grounded import qm, tracefree, word, braid

TAU = 2 * math.pi
PI = math.pi

# --- base tangles in the compass slot order (NW, SW, SE, NE) -------------------------
# 0-tangle: horizontal strands NW-NE (m1) and SW-SE (m2)  -> theta = 0 edge
ZERO = [[(0, 1)], [(1, 1)], [(1, 1)], [(0, 1)]]
# infinity-tangle: vertical strands NW-SW (m1) and NE-SE (m2) -> gamma = 0 edge
INF = [[(0, 1)], [(0, 1)], [(1, 1)], [(1, 1)]]


def east_twists(n, base=None):
    """Integer tangle T_n: n horizontal twists at the East punctures = slots (x3,x4),
    0-based braid index 2. Positive n uses inv=False (derived slope +n)."""
    W = [list(w) for w in (base or ZERO)]
    for _ in range(abs(n)):
        W = braid(W, 2, inv=(n < 0))
    return W


def south_twists(n, base=None):
    """Vertical tangle Q(1/n): n twists at the South punctures = slots (x2,x3),
    0-based braid index 1. Positive n uses inv=True (derived slope +1/n)."""
    W = [list(w) for w in (base or INF)]
    for _ in range(abs(n)):
        W = braid(W, 1, inv=(n > 0))
    return W


def rotate(words):
    """90-degree CCW box rotation: new (NW,SW,SE,NE) = old (NE,NW,SW,SE).
    On curves: (gamma,theta) -> (-theta,gamma); fractions F -> -1/F."""
    return [words[3], words[0], words[1], words[2]]


# --- curve extraction: boundary words -> iota-invariant curve on the double cover T^2 -
def betas(words, eta):
    """Exact beta-angles of the four boundary meridians at parameter eta.
    Generators m1 = i, m2 = cos(eta) i + sin(eta) j; every boundary word of a
    braid-built tangle is a conjugate of a generator, hence a pure quaternion in
    the i-j plane: beta = atan2(j-comp, i-comp). Asserts coplanarity + product = 1."""
    gens = [tracefree((1, 0, 0)), tracefree((math.cos(eta), math.sin(eta), 0))]
    Xs = [word(gens, w) for w in words]
    prod = (1, 0, 0, 0)
    for X in Xs:
        assert abs(X[0]) < 1e-9 and abs(X[3]) < 1e-9, "boundary meridian not coplanar-traceless"
        prod = qm(prod, X)
    assert abs(prod[0] - 1) < 1e-9, "boundary relation x1 x2 x3 x4 = 1 violated"
    return [math.atan2(X[2], X[1]) for X in Xs]


def _wrap(x):
    """wrap to (-pi, pi]"""
    return x - TAU * math.floor(x / TAU + 0.5)


def curve(words, N=360):
    """The tangle's pillowcase curve as polylines on T^2 = [0,2pi)^2: the eta-arc
    (locally unwrapped, continuous coordinates) plus its iota-image."""
    arc = []
    g_prev = t_prev = None
    for k in range(0, N + 1):   # include eta = 0, pi: the corner (reducible) endpoints
        eta = PI * k / N
        b = betas(words, eta)
        g, t = (b[1] - b[0]) % TAU, (b[3] - b[0]) % TAU
        if g_prev is not None:  # continue the local lift
            g = g_prev + _wrap(g - g_prev % TAU)
            t = t_prev + _wrap(t - t_prev % TAU)
        arc.append((g, t))
        g_prev, t_prev = g, t
    iarc = [(-g, -t) for (g, t) in arc]
    return [arc, iarc]


# --- segments + the Conway-sum fiber product ------------------------------------------
def segments(arcs, eps=1e-12):
    """Polylines -> gamma-monotone segments (g0,g1,t0,t1), g0<g1, g0 in [0,2pi),
    split at the gamma-seam so overlap logic is interval arithmetic on [0,2pi]."""
    segs, vertical = [], 0
    for arc in arcs:
        for (p, q) in zip(arc, arc[1:]):
            g0, t0 = p
            g1, t1 = q
            if abs(g1 - g0) < eps:
                vertical += 1
                continue
            if g1 < g0:
                g0, g1, t0, t1 = g1, g0, t1, t0
            shift = TAU * math.floor(g0 / TAU)
            g0, g1 = g0 - shift, g1 - shift
            if g1 <= TAU + eps:
                segs.append((g0, min(g1, TAU), t0, t1))
            else:  # split at the seam
                tm = t0 + (t1 - t0) * (TAU - g0) / (g1 - g0)
                segs.append((g0, TAU, t0, tm))
                segs.append((0.0, g1 - TAU, tm, t1))
    return segs, vertical


def fiber_sum(segsA, segsB, eps=1e-12):
    """The Conway sum on curves: fiber product over gamma, theta adds (sec 23)."""
    out = []
    B = sorted(segsB)
    for (a0, a1, s0, s1) in segsA:
        for (b0, b1, u0, u1) in B:
            if b0 >= a1:
                break
            g0, g1 = max(a0, b0), min(a1, b1)
            if g1 - g0 < eps:
                continue
            def tA(g): return s0 + (s1 - s0) * (g - a0) / (a1 - a0)
            def tB(g): return u0 + (u1 - u0) * (g - b0) / (b1 - b0)
            out.append((g0, g1, tA(g0) + tB(g0), tA(g1) + tB(g1)))
    return out


def tangle_sum(wordsA_or_segs, wordsB_or_segs, N=360):
    """Convenience: accept boundary words or precomputed segments."""
    def to_segs(x):
        return x if (x and isinstance(x[0][0], float)) else segments(curve(x, N))[0]
    return fiber_sum(to_segs(wordsA_or_segs), to_segs(wordsB_or_segs))


# --- seam fiber circles (Smith Thm 4.22, derived independently in sec 24) -------------
def seam_points(segs, tol=1e-9, dedup=1e-6):
    """Seam-touching theta-values of a curve on T^2, folded to P-coordinates:
    {0.0: [...], pi: [...]}, theta in the OPEN interval (0,pi) -- corner points are
    excluded (reducible for rational summands: the S^1 of gluings is all-conjugate,
    no circle; derivation in RESEARCH_LOG sec 24)."""
    out = {0.0: [], PI: []}
    for (g0, g1, t0, t1) in segs:
        for s in (0.0, PI, TAU):
            if g0 - tol <= s <= g1 + tol:
                th = t0 if g1 <= g0 else t0 + (t1 - t0) * (min(max(s, g0), g1) - g0) / (g1 - g0)
                key = PI if s == PI else 0.0
                thP = th % TAU
                if thP > PI:
                    thP = TAU - thP
                if thP < dedup or thP > PI - dedup:
                    continue  # corner
                if all(abs(thP - x) > dedup for x in out[key]):
                    out[key].append(thP)
    for k in out:
        out[k].sort()
    return out


def fiber_circles(segsA, segsB):
    """The seam fiber circles of the Conway sum (the components the generic fiber
    product misses): ONE circle per seam s in {0,pi} per pair of non-corner seam
    points (s,th1) in L1, (s,th2) in L2.  Derivation (sec 24): at the seam the glued
    axis pair is aligned, its stabilizer is S^1, and for irreducible rho_i the S^1 of
    gluings gives distinct characters. Outer North angle along the family:
        cos th_out(alpha) = cos th1 cos th2 + sin th1 sin th2 cos 2alpha,
    so the circle maps 2:1 onto the seam arc [ |th1-th2|, pi-|pi-(th1+th2)| ], its
    fold endpoints being exactly the two iota-lift fiber-product points (coplanar
    gluings). Each circle is returned with its parametrization data."""
    circles = []
    spA, spB = seam_points(segsA), seam_points(segsB)
    for s in (0.0, PI):
        for t1 in spA[s]:
            for t2 in spB[s]:
                circles.append(dict(seam=s, th1=t1, th2=t2,
                                    lo=abs(t1 - t2), hi=PI - abs(PI - (t1 + t2))))
    return circles


def circle_image(c, N=90):
    """P-image samples of a fiber circle over alpha in [0,pi] (cos 2alpha covers
    [-1,1]; the other half repeats the arc -- the 2:1 fold)."""
    pts = []
    for k in range(N + 1):
        ca = math.cos(c['th1']) * math.cos(c['th2']) + \
             math.sin(c['th1']) * math.sin(c['th2']) * math.cos(2 * PI * k / N)
        pts.append((c['seam'], math.acos(max(-1.0, min(1.0, ca)))))
    return pts


def fiber_sum_full(segsA, segsB):
    """Complete Conway sum on curves: generic fiber product + seam fiber circles."""
    return fiber_sum(segsA, segsB), fiber_circles(segsA, segsB)


# --- independent 3D validation engine: explicit quaternion gluing ---------------------
# No coplanarity assumed anywhere below: representations live in full SU(2), the
# gluing rotation is constructed explicitly, and the outer boundary is read off by
# finding the common plane -- so it also *tests* the sec-23 coplanarity lemma.

def _axis(X):
    n = math.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2)
    return (X[1] / n, X[2] / n, X[3] / n)

def _dot3(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def _cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

def _ang3(a, b):
    c = _cross3(a, b)
    return math.atan2(math.sqrt(_dot3(c, c)), _dot3(a, b))  # in [0,pi]

def _quat(axis, ang):
    s = math.sin(ang / 2)
    return (math.cos(ang / 2), s * axis[0], s * axis[1], s * axis[2])

def _adj(g, X):
    """Ad_g X = g X g^{-1}."""
    w, x, y, z = g
    return qm(qm(g, X), (w, -x, -y, -z))

def _rot_between(a, b):
    """A rotation (unit quaternion) carrying axis a to axis b."""
    c = _cross3(a, b)
    s = math.sqrt(_dot3(c, c))
    d = _dot3(a, b)
    if s < 1e-12:
        if d > 0:
            return (1.0, 0.0, 0.0, 0.0)
        p = (1.0, 0.0, 0.0) if abs(a[0]) < 0.9 else (0.0, 1.0, 0.0)
        n = _cross3(a, p)
        nn = math.sqrt(_dot3(n, n))
        return _quat((n[0] / nn, n[1] / nn, n[2] / nn), PI)
    n = (c[0] / s, c[1] / s, c[2] / s)
    return _quat(n, math.atan2(s, d))

def _align_pair(a2, b2, a1, b1):
    """The rotation carrying the ordered axis pair (a2,b2) to (a1,b1) (equal pair
    angles assumed): align a2->a1, then rotate about a1 to bring b2 onto b1."""
    g1 = _rot_between(a2, a1)
    b2p = _axis(_adj(g1, (0.0,) + b2))
    u = tuple(b2p[i] - _dot3(b2p, a1) * a1[i] for i in range(3))
    v = tuple(b1[i] - _dot3(b1, a1) * a1[i] for i in range(3))
    delta = math.atan2(_dot3(_cross3(u, v), a1), _dot3(u, v))
    return qm(_quat(a1, delta), g1)

def _eval4(words, eta):
    gens = [tracefree((1, 0, 0)), tracefree((math.cos(eta), math.sin(eta), 0))]
    return [word(gens, w) for w in words]

FLIP = (0.0, 1.0, 0.0, 0.0)  # conjugation by i: the iota-lift (beta -> -beta)


def pillowcase_read(axes):
    """(gamma, theta) from four coplanar meridian axes, robust at seams (where
    grounded.pillowcase_coords degenerates). Finds the common plane from the most
    transverse axis pair; ASSERTS coplanarity (the sec-23 lemma, here a test)."""
    best, bn = None, 0.0
    for i in range(4):
        for j in range(i + 1, 4):
            c = _cross3(axes[i], axes[j])
            s = math.sqrt(_dot3(c, c))
            if s > bn:
                bn, best = s, tuple(x / s for x in c)
    if best is None:
        return None  # all four axes parallel: corner
    for a in axes:
        assert abs(_dot3(a, best)) < 1e-6, "outer boundary rep NOT coplanar"
    e1 = tuple(axes[0][i] - _dot3(axes[0], best) * best[i] for i in range(3))
    n1 = math.sqrt(_dot3(e1, e1))
    e1 = tuple(x / n1 for x in e1)
    e2 = _cross3(best, e1)
    b = [math.atan2(_dot3(a, e2), _dot3(a, e1)) for a in axes]
    return ((b[1] - b[0]) % TAU, (b[3] - b[0]) % TAU)


def _psi_west(words):
    def f(eta):
        Y = _eval4(words, eta)
        return _ang3(_axis(Y[0]), _axis(Y[1]))
    return f


def _match_eta2(psi2, table, psi_target):
    """All eta2 in (0,pi) with psi2(eta2) = psi_target, from a precomputed
    (eta, psi) table by sign-change scan + live bisection."""
    roots = []
    for (e0, v0), (e1, v1) in zip(table, table[1:]):
        f0, f1 = v0 - psi_target, v1 - psi_target
        if f0 * f1 > 0:
            continue
        a, b_, fa = e0, e1, f0
        for _ in range(60):
            m = (a + b_) / 2
            fm = psi2(m) - psi_target
            if fa * fm <= 0:
                b_ = m
            else:
                a, fa = m, fm
        r = (a + b_) / 2
        if all(abs(r - x) > 1e-7 for x in roots):
            roots.append(r)
    return roots


def glue3d_generic(words1, words2, N=120, M=2400):
    """Independent Conway sum: explicit SO(3) gluing at generic gamma. For each
    eta1, find eta2 with matching glued-pair angle, build the aligning rotation for
    both iota-lifts of rho2, read the outer boundary in full 3D."""
    psi2 = _psi_west(words2)
    table = [(PI * k / M, psi2(PI * k / M)) for k in range(1, M)]
    pts = []
    for k1 in range(1, N):
        eta1 = PI * k1 / N
        X = _eval4(words1, eta1)
        a1, b1 = _axis(X[3]), _axis(X[2])          # glued pair (NE1, SE1)
        psi1 = _ang3(a1, b1)
        if psi1 < 5e-2 or psi1 > PI - 5e-2:
            continue                                # seams handled by glue3d_seam
        for eta2 in _match_eta2(psi2, table, psi1):
            for flip in (False, True):
                Y = _eval4(words2, eta2)
                if flip:
                    Y = [_adj(FLIP, y) for y in Y]
                g = _align_pair(_axis(Y[0]), _axis(Y[1]), a1, b1)   # (NW2,SW2)->(NE1,SE1)
                assert _ang3(_axis(_adj(g, Y[0])), a1) < 1e-6
                assert _ang3(_axis(_adj(g, Y[1])), b1) < 1e-6
                outer = [_axis(X[0]), _axis(X[1]),
                         _axis(_adj(g, Y[2])), _axis(_adj(g, Y[3]))]
                pt = pillowcase_read(outer)
                if pt:
                    pts.append(pt)
    return pts


def _extremal_etas(psi, M=2400):
    """eta values in (0,pi) where the continuous pair-angle psi(eta) in [0,pi]
    touches 0 or pi (tangential extrema: scan + ternary refine)."""
    out = {0.0: [], PI: []}
    vals = [(PI * k / M, psi(PI * k / M)) for k in range(1, M)]
    for target in (0.0, PI):
        lo = (target == 0.0)
        for i in range(1, len(vals) - 1):
            is_min = vals[i][1] <= vals[i - 1][1] and vals[i][1] <= vals[i + 1][1]
            is_max = vals[i][1] >= vals[i - 1][1] and vals[i][1] >= vals[i + 1][1]
            if not (is_min if lo else is_max):
                continue
            a, b_ = vals[i - 1][0], vals[i + 1][0]
            for _ in range(80):  # ternary refine of the extremum
                m1, m2 = a + (b_ - a) / 3, b_ - (b_ - a) / 3
                if (psi(m1) < psi(m2)) == lo:
                    b_ = m2
                else:
                    a = m1
            r = (a + b_) / 2
            if abs(psi(r) - target) < 1e-5 and all(abs(r - x) > 1e-6 for x in out[target]):
                out[target].append(r)
    return out


def _psi_east(words):
    def f(eta):
        X = _eval4(words, eta)
        return _ang3(_axis(X[3]), _axis(X[2]))
    return f


def glue3d_seam(words1, words2, Nalpha=180):
    """Independent seam-circle computation: at eta1*, eta2* where the glued-pair
    angles hit 0 or pi, sweep the S^1 stabilizer explicitly. Returns, per circle,
    the sampled outer thetas and the triple-product pass-sign sequence."""
    circles = []
    e1s, e2s = _extremal_etas(_psi_east(words1)), _extremal_etas(_psi_west(words2))
    for s in (0.0, PI):
        for eta1 in e1s[s]:
            X = _eval4(words1, eta1)
            n = _axis(X[3])                       # the common glued axis
            th1 = _ang3(_axis(X[0]), n)           # angle(NW1, glued axis)
            for eta2 in e2s[s]:
                Y = _eval4(words2, eta2)
                th2 = _ang3(_axis(Y[3]), _axis(Y[0]))
                g0 = _rot_between(_axis(Y[0]), n)
                gam = _ang3(_axis(X[0]), _axis(X[1]))
                assert abs(gam - s) < 1e-4, "outer West angle not on the seam"
                thetas, signs = [], []
                for k in range(Nalpha):
                    g = qm(_quat(n, TAU * k / Nalpha), g0)
                    oNE = _axis(_adj(g, Y[3]))
                    trip = _dot3(_cross3(_axis(X[0]), oNE), n)  # out-of-plane pass sign
                    thetas.append(_ang3(_axis(X[0]), oNE))
                    signs.append(math.copysign(1.0, trip) if abs(trip) > 1e-9 else 0.0)
                circles.append(dict(seam=s, th1=th1, th2=th2, thetas=thetas, signs=signs))
    return circles


# --- validation helpers ----------------------------------------------------------------
def cloud(segs):
    pts = []
    for (g0, g1, t0, t1) in segs:
        for lam in (0.0, 0.5, 1.0):
            pts.append(((g0 + lam * (g1 - g0)) % TAU, (t0 + lam * (t1 - t0)) % TAU))
    return pts


def subtorus_residual(pts, p, q):
    """Max deviation of the cloud from the primitive (q,p) subtorus line
    {q*theta - p*gamma = 0 mod 2pi} (= the slope-p/q curve through the corners)."""
    r = 0.0
    for (g, t) in pts:
        x = (q * t - p * g) % TAU
        r = max(r, min(x, TAU - x))
    return r


def hausdorff(ptsA, ptsB):
    """Symmetric Hausdorff distance on the flat torus."""
    def d(u, v):
        return math.hypot(min(abs(u[0] - v[0]), TAU - abs(u[0] - v[0])),
                          min(abs(u[1] - v[1]), TAU - abs(u[1] - v[1])))
    def h(A, B):
        return max(min(d(a, b) for b in B) for a in A)
    return max(h(ptsA, ptsB), h(ptsB, ptsA))


def sheets_at(segs, g_probe, tol=1e-6):
    """Number of distinct theta values of the curve over gamma = g_probe."""
    vals = []
    for (g0, g1, t0, t1) in segs:
        if g0 - tol <= g_probe <= g1 + tol and g1 > g0:
            t = (t0 + (t1 - t0) * (g_probe - g0) / (g1 - g0)) % TAU
            if all(min(abs(t - v), TAU - abs(t - v)) > tol for v in vals):
                vals.append(t)
    return len(vals)


# --- the battery -----------------------------------------------------------------------
if __name__ == "__main__":
    TOL = 1e-7
    results = []

    def check(name, val, expect_zero=True, tol=TOL):
        ok = (val < tol) if expect_zero else bool(val)
        results.append(ok)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}  (residual {val:.2e})"
              if expect_zero else f"  [{'PASS' if ok else 'FAIL'}] {name}")
        return ok

    print("== compass-pinned tangle calculus: validation battery (RESEARCH_LOG sec 23) ==\n")

    print("(1) base tangles on their edges")
    check("0-tangle on theta=0 (slope 0/1)", subtorus_residual(cloud(segments(curve(ZERO))[0]), 0, 1))
    infsegs, nvert = segments(curve(INF))
    check("inf-tangle on gamma=0 (all segments vertical)", 0.0 if not infsegs and nvert > 0 else 1.0)

    print("(2) East twists: T_n -> slope n  [fraction n]")
    for n in (1, 2, 3, 4, 5, -1, -3):
        r = subtorus_residual(cloud(segments(curve(east_twists(n)))[0]), n, 1)
        check(f"T_{n:+d}: theta = {n}*gamma", r)

    print("(3) South twists: Q(1/n) -> slope 1/n  [fraction 1/n]")
    for n in (2, 3, 5, -2, -3):
        r = subtorus_residual(cloud(segments(curve(south_twists(n)))[0]), 1, n)
        check(f"Q(1/{n}): {n}*theta = gamma", r)

    print("(4) South twists on the 0-tangle are inert (R1 kinks on one strand)")
    r = subtorus_residual(cloud(segments(curve(south_twists(3, base=ZERO)))[0]), 0, 1)
    check("south^3(0-tangle) still on theta=0", r)

    print("(5) rotation: rotate = cyclic slot shift, F -> -1/F")
    for n in (2, 3):
        r = subtorus_residual(cloud(segments(curve(rotate(east_twists(n))))[0]), 1, -n)
        check(f"rotate(T_{n}) on slope -1/{n}", r)
        hd = hausdorff(cloud(segments(curve(rotate(east_twists(n))))[0]),
                       cloud(segments(curve(south_twists(-n)))[0]))
        check(f"rotate(T_{n}) == Q(-1/{n}) as curves (Hausdorff)", hd, tol=2e-2)

    print("(6) Conway sum: T_m + T_n = T_{m+n}  [the sec-22 failure, now derived]")
    for (m, n) in ((1, 1), (2, 3), (-1, 3), (2, -5), (1, -1)):
        s = tangle_sum(east_twists(m), east_twists(n))
        r = subtorus_residual(cloud(s), m + n, 1)
        check(f"T_{m:+d} + T_{n:+d} on slope {m+n}", r)
        hd = hausdorff(cloud(s), cloud(segments(curve(east_twists(m + n)))[0]))
        check(f"T_{m:+d} + T_{n:+d} == T_{m+n:+d} pointwise (Hausdorff)", hd, tol=2e-2)

    print("(7) mixed sum: T_m + Q(1/n) -> slope (mn+1)/n; det of numerator closure = mn+1")
    for (m, n) in ((1, 2), (2, 3), (3, 2)):
        s = tangle_sum(east_twists(m), south_twists(n))
        r = subtorus_residual(cloud(s), m * n + 1, n)
        check(f"T_{m} + Q(1/{n}) on slope {m*n+1}/{n}  [det {m*n+1}]", r)

    print("(8) the P(-2,3,5) sub-assembly: Q(1/3)+Q(1/5), then +Q(-1/2)")
    s35 = tangle_sum(south_twists(3), south_twists(5))
    check("Q(1/3)+Q(1/5) on slope 8/15", subtorus_residual(cloud(s35), 8, 15))
    check("  15 sheets over generic gamma", 0.0 if sheets_at(s35, 1.234) == 15 else 1.0)
    pretzel = fiber_sum(s35, segments(curve(south_twists(-2)))[0])
    check("P(-2,3,5) = +Q(-1/2) on slope 1/30", subtorus_residual(cloud(pretzel), 1, 30))
    check("  30 sheets over generic gamma", 0.0 if sheets_at(pretzel, 1.234) == 30 else 1.0)
    print("  (slope 1/30: numerator-closure det = 1 = det T(3,5) -- the torus-pretzel")
    print("   coincidence of sec 19; the unperturbed pairing is degenerate, hence Smith")
    print("   perturbs before summing. Next step: shear + earring on THIS assembly.)")

    # ------------------------------------------------------------------ seam circles
    def d_P(u, v):
        def dt(a, b):
            x = abs(a - b) % TAU
            return min(x, TAU - x)
        d1 = math.hypot(dt(u[0], v[0]), dt(u[1], v[1]))
        d2 = math.hypot(dt(u[0], -v[0]), dt(u[1], -v[1]))
        return min(d1, d2)

    def hausdorff_P(A, B):
        def h(A_, B_):
            return max(min(d_P(a, b) for b in B_) for a in A_)
        return max(h(A, B), h(B, A))

    print("\n(9) seam fiber circles: counts and arcs (analytic rule)")
    segQ = {n: segments(curve(south_twists(n)))[0] for n in (2, 3, 5)}
    for (name, sA, sB, expect) in (
            ("T_2 + T_3", segments(curve(east_twists(2)))[0], segments(curve(east_twists(3)))[0], 0),
            ("T_2 + Q(1/3)", segments(curve(east_twists(2)))[0], segQ[3], 0),
            ("Q(1/3)+Q(1/5)", segQ[3], segQ[5], 4),
            ("Q(1/3)+Q(1/3)  [the 8_19 assembly]", segQ[3], segQ[3], 2)):
        circ = fiber_circles(sA, sB)
        check(f"{name}: {len(circ)} circles (expect {expect})",
              0.0 if len(circ) == expect else 1.0)
    circ35 = fiber_circles(segQ[3], segQ[5])
    expected_arcs = sorted([(0.0, 4 * PI / 15, 14 * PI / 15), (0.0, 2 * PI / 15, 8 * PI / 15),
                            (PI, 2 * PI / 15, 8 * PI / 15), (PI, 4 * PI / 15, 14 * PI / 15)])
    got_arcs = sorted([(c['seam'], c['lo'], c['hi']) for c in circ35])
    err = max(max(abs(a - b) for a, b in zip(x, y)) for x, y in zip(expected_arcs, got_arcs))
    check("Q(1/3)+Q(1/5) arcs = {[4pi/15,14pi/15],[2pi/15,8pi/15]} x both seams", err)
    circ33 = fiber_circles(segQ[3], segQ[3])
    err = max(abs(c['lo']) for c in circ33) + max(abs(c['hi'] - 2 * PI / 3) for c in circ33)
    check("Q(1/3)+Q(1/3) arcs = [0, 2pi/3] both seams (corner-touching!)", err)
    # continuity: circle fold-endpoints are seam values of the generic fiber product
    s35 = tangle_sum(south_twists(3), south_twists(5))
    seamvals = seam_points(s35, dedup=1e-4)
    err = 0.0
    for c in circ35:
        for e in (c['lo'], c['hi']):
            err = max(err, min(abs(e - v) for v in seamvals[c['seam']]))
    check("circle endpoints = generic-sum seam values (coplanar-gluing continuity)", err, tol=1e-4)

    print("(10) independent 3D gluing engine vs analytic fiber product (generic gamma)")
    for (name, w1, w2, direct) in (
            ("T_2 + T_3", east_twists(2), east_twists(3), east_twists(5)),):
        pts3d = glue3d_generic(w1, w2, N=90)
        ref = cloud(segments(curve(direct))[0])
        hd = hausdorff_P(pts3d, [p for p in ref if 0.20 < p[0] % PI < PI - 0.05])
        check(f"{name}: 3D-glued cloud == T_5 curve (Hausdorff, away from seams)", hd, tol=0.12)
    pts3d = glue3d_generic(south_twists(3), south_twists(5), N=120)
    ref = [p for p in cloud(s35) if 0.20 < p[0] % PI < PI - 0.05]
    hd = hausdorff_P(pts3d, ref)
    check("Q(1/3)+Q(1/5): 3D-glued cloud == analytic fiber product", hd, tol=0.12)

    print("(11) independent 3D seam sweep vs the analytic circle rule")
    sw = glue3d_seam(south_twists(3), south_twists(5), Nalpha=180)
    check(f"3D finds {len(sw)} circles (expect 4)", 0.0 if len(sw) == 4 else 1.0)
    worst_pair, worst_fn, worst_sc = 0.0, 0.0, 0
    for c3 in sw:
        best = min(circ35, key=lambda c: abs(c['seam'] - c3['seam'])
                   + abs(c['th1'] - c3['th1']) + abs(c['th2'] - c3['th2']))
        worst_pair = max(worst_pair, abs(best['th1'] - c3['th1']) + abs(best['th2'] - c3['th2']))
        c0 = math.cos(c3['th1']) * math.cos(c3['th2'])
        s0 = math.sin(c3['th1']) * math.sin(c3['th2'])
        Na = len(c3['thetas'])
        d0 = 0.0 if abs(c3['thetas'][0] - best['lo']) < abs(c3['thetas'][0] - best['hi']) else PI
        fn = max(abs(math.cos(t) - (c0 + s0 * math.cos(TAU * k / Na + d0)))
                 for k, t in enumerate(c3['thetas']))
        worst_fn = max(worst_fn, fn)
        sgn = [s for s in c3['signs'] if s != 0.0]
        flips = sum(1 for a, b in zip(sgn, sgn[1:] + sgn[:1]) if a != b)
        worst_sc = max(worst_sc, abs(flips - 2))
        worst_sc = max(worst_sc, 0 if min(c3['thetas']) >= best['lo'] - 1e-3
                       and max(c3['thetas']) <= best['hi'] + 1e-3 else 1)
    check("  circle (th1,th2) data matches analytic pairs", worst_pair, tol=1e-5)
    check("  theta_out(alpha) == spherical-cosine formula, all samples", worst_fn, tol=1e-9)
    check("  pass-sign flips exactly twice per loop (2:1 fold), arcs within [lo,hi]",
          float(worst_sc), tol=0.5)
    sw33 = glue3d_seam(south_twists(3), south_twists(3), Nalpha=180)
    check(f"8_19 assembly: 3D finds {len(sw33)} circles (expect 2), corner-touching",
          0.0 if len(sw33) == 2 and max(min(c['thetas']) for c in sw33) < 0.08 else 1.0)

    print(f"\n{'ALL PASS' if all(results) else 'FAILURES PRESENT'} "
          f"({sum(results)}/{len(results)})")
