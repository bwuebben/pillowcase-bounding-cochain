r"""Floer complexes of the sheared curves of Section 4, computed numerically (Appendix A.1(v), summarized in Remark 4.21).

The perturbed variety W_eps of the tangle of Hedden, Herald and Kirk for T(3,5), (r, s) = (2, -1), is traced with the
continuation code of the companion paper on T(3,n) (module t3n_pillowcase).  Everything downstream is computed here,
independently of that module's complex code:

  * the restriction map to the pillowcase, from the quaternion representation [HHK2, (31)], and a continuous lift
    to R^2, compared with the lift returned by the continuation code;
  * the shear S_{-m}(gamma, theta) = (gamma, theta + (q-5) gamma) of Lemma 4.6;
  * the lifts S^k_+- of the earring curve L_0 (Section 2.2) as graphs over gamma, and the generators as the
    crossings of the sheared curve with these lifts;
  * degrees from (2.3): deg x^+ = 2k + mu and deg x^- = deg x^+ - 1, where mu is the Maslov index of the sub-path from
    r_+, computed on the unsheared curve against the line field spanned by (1, 6-q).  By Lemma 4.6(f) this equals the
    index of the sheared sub-path against the line field of slope one.  (An index obtained by unwrapping tangent
    angles along the sheared polygonal curve is unreliable near the junctions when eps is small or q is large,
    because a single step can turn by more than pi/2 there.)
  * bigons: for two generators on the same lift of L_0, the loop formed by the arc of that lift and the sub-path of the curve is
    accepted if it is a simple closed curve, encloses no lattice point, is oriented counterclockwise and has convex
    corners (that is, if its lift bounds an embedded disc in R^2 \ (pi Z)^2);
  * graded homology over F_2.

Floating-point numerics; these computations are not used in the proofs.
"""
import math

import numpy as np

PI, TAU = math.pi, 2 * math.pi


# ------------------------------------------------------------------ restriction map and lift
def qmul(p, q):
    w1, x1, y1, z1 = p.T
    w2, x2, y2, z2 = q.T
    return np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], axis=1)


def qexp(ang, ax):
    return np.concatenate([np.cos(ang)[:, None], np.sin(ang)[:, None] * ax], axis=1)


def boundary(X, p, q, r, s, eA, eB):
    """The boundary representation (a, b, c) at points X = (u, v, tau) of W_eps [HHK2, (31)-(32)]."""
    u, v, tau = X[:, 0], X[:, 1], np.clip(X[:, 2], -1, 1)
    QA = np.tile([1.0, 0, 0], (len(u), 1))
    QB = np.stack([tau, np.sqrt(np.maximum(0, 1 - tau ** 2)), 0 * tau], axis=1)
    A1 = lambda e: qexp(e * u, QA)
    A2 = lambda e: qexp(e * eA * np.sin(u), QA)
    B1 = lambda e: qexp(e * v, QB)
    B2 = lambda e: qexp(e * eB * np.sin(v), QB)
    a = qmul(qmul(A1(s + p), A2(q - r)), qmul(B1(q - r), B2(-(s + p))))
    b = qmul(qmul(B1(-r), B2(-s)), qmul(A1(s), A2(-r)))
    x = qmul(B2(s), B1(r))
    c = qmul(qmul(x * np.array([1, -1, -1, -1]), a), x)
    return a, b, c


def pillow_from(a, b, c):
    """The point (gamma, theta) with a -> i, b -> e^{gamma k} i, c -> e^{theta k} i up to conjugation."""
    av, bv, cv = a[:, 1:], b[:, 1:], c[:, 1:]
    cg = np.clip((av * bv).sum(1), -1, 1)
    g = np.arccos(cg)
    sg = np.sin(g)
    e2 = (bv - cg[:, None] * av) / np.where(sg > 1e-12, sg, 1)[:, None]
    cp = cv - (cv * av).sum(1)[:, None] * av
    e2 = np.where((sg > 1e-12)[:, None], e2, cp / np.maximum(np.linalg.norm(cp, axis=1), 1e-15)[:, None])
    th = np.arctan2((cv * e2).sum(1), (cv * av).sum(1))
    return np.c_[g, th]


def pillow(X, p, q, r, s, eA, eB):
    a, b, c = boundary(X, p, q, r, s, eA, eB)
    assert np.abs(a[:, 0]).max() < 1e-7 and np.abs(b[:, 0]).max() < 1e-7, 'not traceless'
    return pillow_from(a, b, c)


def lift(P):
    """A continuous lift to R^2 of a sampled path in the pillowcase R^2/G."""
    out = np.empty_like(P)
    out[0] = P[0]
    pg, pt = P[0]
    for i in range(1, len(P)):
        best = None
        for sg in (1.0, -1.0):
            G0, T0 = sg * P[i, 0], sg * P[i, 1]
            cg = G0 + TAU * round((pg - G0) / TAU)
            ct = T0 + TAU * round((pt - T0) / TAU)
            d = (cg - pg) ** 2 + (ct - pt) ** 2
            if best is None or d < best[0]:
                best = (d, cg, ct)
        pg, pt = best[1], best[2]
        out[i] = (pg, pt)
    return out


def normalize_arc(L):
    """Translate the lift of the arc so that it starts at (0, 0) and runs into gamma > 0."""
    c0 = np.round(L[0] / PI).astype(int)
    if c0[0] % 2 or c0[1] % 2:
        L = L[::-1].copy()
        c0 = np.round(L[0] / PI).astype(int)
    assert c0[0] % 2 == 0 and c0[1] % 2 == 0, 'the arc does not start at a lift of (0,0)'
    L = L - c0 * PI
    if L[min(50, len(L) - 1), 0] < 0:
        L = -L
    L[0] = 0.0
    L[-1] = np.round(L[-1] / PI) * PI
    return L


def shear(L, q):
    """S_{-m}, m = (q-5)/2: (gamma, theta) -> (gamma, theta + (q-5) gamma)."""
    L = np.array(L, float)
    L[:, 1] = L[:, 1] + (q - 5) * L[:, 0]
    return L


# ------------------------------------------------------------------ the lifts of the earring curve and the generators
def strand_theta(sig, gam, eps):
    """theta of the lift S^0_sig of L_0 over gamma (Section 2.2): S_+ is t -> (t + eps sin t + pi/2, t - eps sin t + pi/2)
    and S_- its image under (gamma, theta) -> (-gamma, -theta)."""
    target = (gam - PI / 2) if sig > 0 else (-gam - PI / 2)
    t = np.array(target, float)
    for _ in range(60):
        t = t - (t + eps * np.sin(t) - target) / (1 + eps * np.cos(t))
    return (gam - 2 * eps * np.sin(t)) if sig > 0 else (gam + 2 * eps * np.sin(t))


def generators(L, eps):
    """Crossings of the polyline L with the lifts S^k_sig of L_0: dicts with sig, k, pos (fractional index), pt, dir."""
    out = []
    for sig in (1, -1):
        f = L[:, 1] - strand_theta(sig, L[:, 0], eps)
        kk = np.floor(f / TAU)
        for i in np.nonzero(kk[:-1] != kk[1:])[0]:
            lo, hi = sorted((f[i], f[i + 1]))
            for k in range(math.ceil(lo / TAU), math.floor(hi / TAU) + 1):
                if not (lo < TAU * k < hi):
                    continue
                a_, b_ = 0.0, 1.0
                fa = f[i] - TAU * k
                for _ in range(60):
                    mid = 0.5 * (a_ + b_)
                    Pm = L[i] + mid * (L[i + 1] - L[i])
                    fm = Pm[1] - strand_theta(sig, Pm[0], eps) - TAU * k
                    if (fm > 0) == (fa > 0):
                        a_, fa = mid, fm
                    else:
                        b_ = mid
                mid = 0.5 * (a_ + b_)
                out.append(dict(sig=sig, k=k, pos=i + mid, pt=L[i] + mid * (L[i + 1] - L[i]),
                                dir=1 if f[i + 1] > f[i] else -1))
    out.sort(key=lambda d: d['pos'])
    return out


# ------------------------------------------------------------------ Maslov index, sub-paths, embeddedness
def mu_passages(L, ref):
    """Signed number of passages of the tangent line of the polyline L through the constant line field at angle ref
    (counterclockwise passages count +1)."""
    d = np.diff(L, axis=0)
    d = d[np.hypot(d[:, 0], d[:, 1]) > 1e-14]
    if len(d) < 2:
        return 0
    a = 2 * np.arctan2(d[:, 1], d[:, 0])
    jumps = np.abs(np.mod(np.diff(a) + PI, TAU) - PI)
    assert jumps.max() < 0.9 * PI, ('unsafe angle unwrap', jumps.max())
    ang = np.unwrap(a) / 2
    f = np.floor((ang - ref - 1e-12) / PI)
    return int(f[-1] - f[0])


def self_crossings(L, cell=0.05):
    """Proper crossings of non-adjacent segments of the polyline L."""
    P, Q = L[:-1], L[1:]
    lo = np.floor(np.minimum(P, Q) / cell).astype(int)
    hi = np.floor(np.maximum(P, Q) / cell).astype(int)
    grid = {}
    for j in range(len(P)):
        for cx in range(lo[j, 0], hi[j, 0] + 1):
            for cy in range(lo[j, 1], hi[j, 1] + 1):
                grid.setdefault((cx, cy), []).append(j)
    found = set()
    for segs in grid.values():
        if len(segs) < 2:
            continue
        S = np.array(segs)
        for ii in range(len(S)):
            i = S[ii]
            js = S[ii + 1:]
            js = js[np.abs(js - i) > 1]
            if len(js) == 0:
                continue
            a, b = P[i], Q[i]
            c, dd = P[js], Q[js]
            r = b - a
            sv = dd - c
            den = r[0] * sv[:, 1] - r[1] * sv[:, 0]
            ok = np.abs(den) > 1e-18
            qp = c - a
            tt = np.where(ok, (qp[:, 0] * sv[:, 1] - qp[:, 1] * sv[:, 0]) / np.where(ok, den, 1), -1)
            uu = np.where(ok, (qp[:, 0] * r[1] - qp[:, 1] * r[0]) / np.where(ok, den, 1), -1)
            hit = ok & (tt > 1e-12) & (tt < 1 - 1e-12) & (uu > 1e-12) & (uu < 1 - 1e-12)
            for j in js[hit]:
                found.add((min(i, j), max(i, j)))
    return sorted(found)


def sub_polyline(L, p0, p1):
    def pt(pos):
        i = min(int(math.floor(pos)), len(L) - 2)
        return L[i] + (pos - i) * (L[i + 1] - L[i])
    if p0 <= p1:
        mid = L[int(math.floor(p0)) + 1:int(math.floor(p1)) + 1]
    else:
        mid = L[int(math.floor(p1)) + 1:int(math.floor(p0)) + 1][::-1]
    return np.vstack([pt(p0), mid, pt(p1)])


def lattice_winding(loop):
    """Nonzero winding numbers of the closed polyline loop around the lattice points (pi Z)^2."""
    lo = np.floor(loop.min(0) / PI).astype(int) - 1
    hi = np.ceil(loop.max(0) / PI).astype(int) + 1
    W = {}
    cl = np.vstack([loop, loop[:1]])
    for I in range(lo[0], hi[0] + 1):
        for J in range(lo[1], hi[1] + 1):
            e = cl - np.array([I * PI, J * PI])
            a = np.arctan2(e[:, 1], e[:, 0])
            w = np.sum(np.mod(np.diff(a) + PI, TAU) - PI) / TAU
            if abs(w) > 1e-3:
                W[(I, J)] = int(round(w))
    return W


def area(loop):
    x, y = loop[:, 0], loop[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def bigon_test(L, G, p, q, eps):
    """Is there an embedded bigon from the generator p to the generator q (same lift of L_0)?  Returns (bigon, blocked,
    winding, area)."""
    g0, g1 = p['pt'][0], q['pt'][0]
    lo_pos, hi_pos = sorted((p['pos'], q['pos']))
    glo, ghi = sorted((g0, g1))
    blocked = any(lo_pos < h['pos'] < hi_pos and glo < h['pt'][0] < ghi
                  for h in G if h['sig'] == p['sig'] and h['k'] == p['k'] and h is not p and h is not q)
    gs = np.linspace(g0, g1, 400)
    strand = np.c_[gs, strand_theta(p['sig'], gs, eps) + TAU * p['k']]
    strand[0], strand[-1] = p['pt'], q['pt']
    arc = sub_polyline(L, q['pos'], p['pos'])
    loop = np.vstack([strand, arc[1:-1]])
    W = lattice_winding(loop)
    A = area(loop)
    unit = lambda v: v / np.hypot(*v)
    sin_dir = unit(strand[-1] - strand[-2])
    aout = unit(arc[1] - arc[0])
    ain = unit(arc[-1] - arc[-2])
    sout = unit(strand[1] - strand[0])
    convex = (sin_dir[0] * aout[1] - sin_dir[1] * aout[0] > 0) and (ain[0] * sout[1] - ain[1] * sout[0] > 0)
    return (not blocked) and (not W) and A > 0 and convex, blocked, W, A


# ------------------------------------------------------------------ homology
def rank_f2(M):
    M = (np.array(M, dtype=np.int64) % 2).copy()
    r = 0
    for c in range(M.shape[1]):
        piv = next((k for k in range(r, M.shape[0]) if M[k, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for k in range(M.shape[0]):
            if k != r and M[k, c]:
                M[k] ^= M[r]
        r += 1
    return r


def graded_homology(deg, edges, n):
    """Graded F_2 homology; edges (i, j) are bigons from generator i to generator j, deg j = deg i - 1."""
    D = np.zeros((n, n), int)
    for (i, j) in edges:
        D[j, i] ^= 1
    H = []
    for g in range(4):
        idx = [i for i in range(n) if deg[i] == g]
        idx_in = [i for i in range(n) if deg[i] == (g + 1) % 4]
        idx_out = [i for i in range(n) if deg[i] == (g - 1) % 4]
        r_out = rank_f2(D[np.ix_(idx_out, idx)]) if idx and idx_out else 0
        r_in = rank_f2(D[np.ix_(idx, idx_in)]) if idx and idx_in else 0
        H.append(len(idx) - r_out - r_in)
    return tuple(H), rank_f2(D)


def predicted(q):
    """Theorem 4.18(3): (1 + N_3, N_1, N_1, N_3)."""
    N1, N3 = -(-(q + 1) // 4), (q + 1) // 4
    return (1 + N3, N1, N1, N3)


def Nq(q):
    return sum(1 for j in range(0, 200) if (q - 6) / 12 < j < 5 * (q - 6) / 12)
