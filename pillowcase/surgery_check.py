#!/usr/bin/env python3
"""
surgery_check.py -- the combinatorial surgery lemma, machine-checked
(correspondence note notes/cochain-search-note, Part II "Reduction 1"; RESUME sec 0d).

CLAIM (combinatorial surgery lemma, deficit members). Let b be the computed
bounding cochain, supported on P self-crossings of blue. Smooth blue at each
support crossing -- iota-equivariantly on the T^2 cover (both preimages), in one
of the two sector pairs -- to get a curve (or curves) blue'. Then for the correct
sector choice,

    bigon_matrix(red, blue')  ==  partial_b(red, blue)    entrywise over F_2,

i.e. the b-DEFORMED differential of the immersed pair equals the UNDEFORMED
(bigon-only) differential of the SURGERED pair. Mechanism: a triangle through s
becomes a bigon of blue' rounding the smoothed corner; the q=5 quadrilateral
through {s_A,s_B} becomes a bigon passing both smoothing sites; a bigon of d
cancelled by the deformation reappears PAIRED with a polygon-turned-bigon and
dies mod 2.

Per member (q=5 with the build_geometry perturbation, matching the certified
b2_result numbers; q=7 with the generic family perturbation):
  1. reproduce the naive complex (gens, d);
  2. recompute partial_b independently (targeted triangle sweep via
     deform.triangle_contributions only_s, P-union over preimages; targeted quad
     count via polygons.polygon_through, both preimage pairs, both cyclic
     orders) and gate it against I^natural;
  3. for every iota-consistent sector combination (2 per support crossing):
     smooth, rebuild, recount bigons, compare to partial_b entrywise;
  4. PASS iff at least one combination reproduces partial_b exactly with the
     generator set unchanged point-by-point.
The other sector choice corresponds to the other CF(L,L) generator at the
crossing; its outcome is reported, not gated.

Pure stdlib. Runtime a few minutes. Exit 0 iff ALL PASS.
"""
import sys, math, itertools
from bigons import (intersections_detailed, arc_of, is_lune,
                    contains_vertex, edges_of, _tdist, CORNERS)
from polygons import self_intersections_detailed, polygon_through
from deform import (build_geometry_p, build_pretzel, bigon_matrix,
                    triangle_contributions, rank_f2)
from earring import P_point
from tangles import TAU

FAIL = 0
EPS = 0.006          # geometric radius of each smoothing ball


def check(cond, msg):
    global FAIL
    print(("  [PASS] " if cond else "  [FAIL] ") + msg)
    if not cond:
        FAIL = 1


def iota(p):
    return ((-p[0]) % TAU, (-p[1]) % TAU)


def entries(M):
    n = len(M)
    return sorted((i, j) for i in range(n) for j in range(n) if M[i][j])


# ---------------------------------------------------------------------------
# smoothing machinery
# ---------------------------------------------------------------------------
def pt_on_edge(E, k, t):
    _, a, b = E[k]
    return ((a[0] + t * (b[0] - a[0])) % TAU, (a[1] + t * (b[1] - a[1])) % TAU)


def cut_endpoints(E, c, eps):
    """Four cut points of crossing c at geometric radius ~eps: on branch A
    (edge kA, param tA -/+ d) and branch B. Returns name -> (k, t, pt)."""
    out = {}
    for br, kk, tt in (("A", c['kA'], c['tA']), ("B", c['kB'], c['tB'])):
        _, a, b = E[kk]
        L = math.hypot(b[0] - a[0], b[1] - a[1])
        d = min(eps / L, tt / 2, (1 - tt) / 2)
        assert d > 1e-6, "cut too close to a polyline vertex; reduce EPS"
        for sgn, name in ((-1, br + "-"), (+1, br + "+")):
            t = tt + sgn * d
            out[name] = (kk, t, pt_on_edge(E, kk, t))
    return out


def chords_for(ep, pairing):
    """The two smoothing chords. pairing 1: A- <-> B+, B- <-> A+ (rounds one
    sector pair); pairing 2: A- <-> B-, A+ <-> B+ (rounds the other)."""
    if pairing == 1:
        return [("A-", "B+"), ("B-", "A+")]
    return [("A-", "B-"), ("A+", "B+")]


def iota_pairing(E, c1, p1, c2, eps):
    """The pairing at preimage c2 whose chord set is the iota-image of the
    chords of (c1, p1); matched by chord midpoints. Returns (p2, mismatch)."""
    ep1 = cut_endpoints(E, c1, eps)
    targets = []
    for u, v in chords_for(ep1, p1):
        m = ((ep1[u][2][0] + ep1[v][2][0]) / 2, (ep1[u][2][1] + ep1[v][2][1]) / 2)
        targets.append(iota(m))
    ep2 = cut_endpoints(E, c2, eps)
    best = None
    for p2 in (1, 2):
        worst = 0.0
        for u, v in chords_for(ep2, p2):
            m = ((ep2[u][2][0] + ep2[v][2][0]) / 2, (ep2[u][2][1] + ep2[v][2][1]) / 2)
            worst = max(worst, min(_tdist(m, t) for t in targets))
        if best is None or worst < best[1]:
            best = (p2, worst)
    return best


def smooth_curve(poly, cuts, eps):
    """cuts: list of (crossing dict, pairing). Returns the smoothed curve as a
    list of closed polylines (components).

    All cut endpoints are cyclically sorted along the curve by (edge, param);
    the interval between the two endpoints on one branch of one crossing is the
    REMOVED stub, every other interval is a KEPT arc; smoothing chords rejoin
    the endpoints; walking the arc+chord graph yields the components."""
    E = edges_of(poly)
    n_edges = len(poly) - 1
    nodes, chords, removed = [], [], set()
    for c, pairing in cuts:
        ep = cut_endpoints(E, c, eps)
        ids = {}
        for name, (k, t, p) in ep.items():
            ids[name] = len(nodes)
            nodes.append((k, t, p))
        for br in ("A", "B"):
            removed.add(frozenset((ids[br + "-"], ids[br + "+"])))
        for u, v in chords_for(ep, pairing):
            chords.append((ids[u], ids[v]))

    order = sorted(range(len(nodes)), key=lambda i: (nodes[i][0], nodes[i][1]))
    arcs = {}
    for pos in range(len(order)):
        i, j = order[pos], order[(pos + 1) % len(order)]
        if frozenset((i, j)) in removed:
            continue
        ki, ti, _ = nodes[i]
        kj, tj, _ = nodes[j]
        mids = []
        if not (ki == kj and tj >= ti):        # same-edge forward arc: no vertices
            v = (ki + 1) % n_edges
            while True:
                mids.append(poly[v])
                if v == kj:
                    break
                v = (v + 1) % n_edges
        arcs[i] = (j, mids)
        arcs[j] = (i, list(reversed(mids)))
    assert set(arcs) == set(range(len(nodes))), "every cut point needs one kept arc"

    chord_of = {}
    for u, v in chords:
        chord_of[u] = v
        chord_of[v] = u
    assert set(chord_of) == set(range(len(nodes))), "every cut point needs one chord"

    comps, seen = [], set()
    for start in range(len(nodes)):
        if start in seen:
            continue
        pts = [nodes[start][2]]
        cur, via_arc = start, True
        while True:
            seen.add(cur)
            if via_arc:
                nxt, mids = arcs[cur]
                pts.extend(mids)
            else:
                nxt = chord_of[cur]
            pts.append(nodes[nxt][2])
            cur, via_arc = nxt, not via_arc
            if cur == start and via_arc:
                break
        if _tdist(pts[0], pts[-1]) > 1e-9:
            pts.append(pts[0])
        comps.append(pts)
    return comps


# ---------------------------------------------------------------------------
# bigons against a multi-component blue
# ---------------------------------------------------------------------------
def bigon_matrix_multi(red, comps):
    gens = []
    for ci, comp in enumerate(comps):
        for h in intersections_detailed(red, comp):
            h['comp'] = ci
            gens.append(h)
    n = len(gens)
    d = [[0] * n for _ in range(n)]
    ncvs = [[k for k, p in enumerate(comp[:-1])
             if any(_tdist(p, c) < 1e-6 for c in CORNERS)] for comp in comps]
    for i in range(n):
        for j in range(n):
            if i == j or gens[i]['comp'] != gens[j]['comp']:
                continue
            x, y = gens[i], gens[j]
            comp = comps[x['comp']]
            for fa in (True, False):
                alpha = arc_of(red, x['kA'], x['tA'], y['kA'], y['tA'], fa)
                for fb in (True, False):
                    if contains_vertex(comp, y['kB'], x['kB'], fb, ncvs[x['comp']]):
                        continue
                    beta = arc_of(comp, y['kB'], y['tB'], x['kB'], x['tB'], fb)
                    if is_lune(alpha, beta, x['pt'], y['pt']):
                        d[i][j] ^= 1
    return gens, d


# ---------------------------------------------------------------------------
# the deformed differential, recomputed independently
# ---------------------------------------------------------------------------
def deformed_matrix(red, blue, gens, supports):
    """partial_b for b = sum(supports): d + Tri (P-union over preimages per
    support) + (if two supports) Quad over both preimage pairs and both cyclic
    orders (deform_full.build_quad convention)."""
    _, d = bigon_matrix(red, blue)
    n = len(gens)
    M = [row[:] for row in d]
    reps, Tri = triangle_contributions(red, blue, gens,
                                       maxspan_blue=230, maxspan_red=60,
                                       only_s=[pp for pp, _ in supports])
    for si in Tri:
        Ms = Tri[si][1]
        for i in range(n):
            for j in range(n):
                M[i][j] ^= Ms[i][j]
    if len(supports) == 2:
        (ppA, preA), (ppB, preB) = supports
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                cnt = 0
                for sa in preA:
                    for sb in preB:
                        cnt += polygon_through(red, blue, gens[i], gens[j],
                                               [sa, sb], maxspan=230)
                        cnt += polygon_through(red, blue, gens[i], gens[j],
                                               [sb, sa], maxspan=230)
                M[i][j] ^= cnt % 2
    return d, M


# ---------------------------------------------------------------------------
# support location: nearest P-orbit of self-crossings to a target point
# ---------------------------------------------------------------------------
def locate_supports(blue, targets):
    scross = self_intersections_detailed(blue)
    orbits = []                                  # (P_point, [preimages])
    for s in scross:
        pp = P_point(s['pt'])
        for o in orbits:
            if _tdist(o[0], pp) < 3e-3:
                o[1].append(s)
                break
        else:
            orbits.append((pp, [s]))
    out = []
    for tgt in targets:
        cands = sorted((o for o in orbits if _tdist(o[0], tgt) < 0.03),
                       key=lambda o: _tdist(o[0], tgt))
        out.append(cands)
    return out


# ---------------------------------------------------------------------------
# the check, per member
# ---------------------------------------------------------------------------
def run_member(q, builder, support_targets, inat, expect_naive_entries=None):
    print(f"\n== q={q}: surger blue at {len(support_targets)} crossing(s), "
          f"target HF = I^natural = {inat} ==")
    red, blue, _ = builder()
    E = edges_of(blue)
    gens, d0 = bigon_matrix(red, blue)
    n = len(gens)
    naive = n - 2 * rank_f2(d0)
    print(f"naive: {n} gens, bigons {entries(d0)}, HF = {naive}")
    if expect_naive_entries is not None:
        check(entries(d0) == sorted(expect_naive_entries),
              f"naive bigon entries match the certified ones {sorted(expect_naive_entries)}")

    # locate the support: among self-crossing P-orbits near the published
    # coordinates (they drift with the perturbation), the support is the
    # combination whose deformed differential reaches I^natural -- unique at
    # deficit members per the certified search (pretzel_solve).
    cand_lists = locate_supports(blue, support_targets)
    for cl, tgt in zip(cand_lists, support_targets):
        check(len(cl) >= 1 and all(len(o[1]) == 2 for o in cl),
              f"candidate orbits near P{tgt}: "
              f"{[tuple(round(v, 3) for v in o[0]) for o in cl]}, 2 preimages each")
    supports, Mb = None, None
    for combo in itertools.product(*cand_lists):
        if len({id(o) for o in combo}) < len(combo):
            continue
        d, M = deformed_matrix(red, blue, gens, list(combo))
        if n - 2 * rank_f2(M) == inat:
            supports, Mb = list(combo), M
            break
    check(supports is not None,
          "a candidate support reaches deformed HF = I^natural (the certified cochain)")
    if supports is None:
        return []
    print(f"support: {[tuple(round(v, 3) for v in pp) for pp, _ in supports]}; "
          f"partial_b entries {entries(Mb)}, HF = {inat}")

    # clearances: smoothing balls clear of generators, red, and each other
    for pp, pre in supports:
        for c in pre:
            dg = min(_tdist(c['pt'], g['pt']) for g in gens)
            dr = min(_tdist(c['pt'], p) for p in red)
            check(dg > 4 * EPS and dr > 4 * EPS,
                  f"cut at P{tuple(round(v, 3) for v in pp)}: clear of gens "
                  f"({dg:.3f}) and red ({dr:.3f}); 4*EPS = {4 * EPS:.3f}")
    allpre = [c for _, pre in supports for c in pre]
    sep = min(_tdist(a['pt'], b['pt'])
              for i, a in enumerate(allpre) for b in allpre[i + 1:])
    check(sep > 4 * EPS, f"cut sites mutually separated ({sep:.3f} > {4 * EPS:.3f})")

    ok = []
    for combo in itertools.product((1, 2), repeat=len(supports)):
        cuts = []
        for (pp, pre), p1 in zip(supports, combo):
            c1, c2 = pre
            p2, mism = iota_pairing(E, c1, p1, c2, EPS)
            assert mism < 3 * EPS, f"iota-matching failed (mismatch {mism:.4f})"
            cuts.append((c1, p1))
            cuts.append((c2, p2))
        comps = smooth_curve(blue, cuts, EPS)
        gens2, d2 = bigon_matrix_multi(red, comps)
        perm = []
        for g in gens:
            m = [i for i, h in enumerate(gens2) if _tdist(h['pt'], g['pt']) < 1e-3]
            perm.append(m[0] if len(m) == 1 else None)
        same = (len(gens2) == n and all(p is not None for p in perm))
        if same:
            d2p = [[d2[perm[i]][perm[j]] for j in range(n)] for i in range(n)]
            ent, hf = entries(d2p), n - 2 * rank_f2(d2p)
            match = (ent == entries(Mb))
        else:
            ent, hf, match = None, None, False
        if match:
            ok.append(combo)
        print(f"  sectors {combo}: {len(comps)} comp(s), {len(gens2)} gens, "
              f"bigons {ent}, HF = {hf}" + ("   <-- matches partial_b" if match else ""))
    check(len(ok) >= 1,
          "some sector choice reproduces partial_b exactly "
          f"(bigons(red, blue') == partial_b entrywise, HF = {inat})")
    return ok


if __name__ == "__main__":
    print("== surgery_check: partial_b(red, blue) =?= bigons(red, surgered blue) ==")
    # q=5: b = s_A + s_B, quadrilateral mechanism (paper2 Computation 1.3(i));
    # build_geometry perturbation, matching the certified b2_result numbers.
    run_member(5, build_geometry_p, [(0.028, 1.272), (3.057, 4.981)], inat=7,
               expect_naive_entries=[(1, 0), (4, 6)])
    # q=7: b = single crossing, triangle mechanism (Computation 1.3(ii));
    # generic family perturbation red=(0.16, 0.40).
    run_member(7, lambda: build_pretzel(3), [(0.05, 5.41)], inat=9)
    print("\n" + ("ALL PASS -- the combinatorial surgery lemma holds at both "
                  "deficit members" if not FAIL else "FAILURES above"))
    sys.exit(FAIL)
