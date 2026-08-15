#!/usr/bin/env python3
"""
solve_b2.py -- step (e), part 5: finite support search for q=5
(RESEARCH_LOG sec 29/30).

Loads the cached deformed-differential tables (deform_full.json: d, Tri per
crossing, Quad per distinct crossing-pair) and evaluates, for any support
B subset {crossings},

    D_tab[i][j] = d[i][j]  XOR  (+)_{S in B} Tri[S][i][j]
                          XOR  (+)_{ {S,S'} subset B } Quad[{S,S'}][i][j]

then searches all subsets of the ACTIVE crossings (those appearing in Tri or Quad)
for rank(D_tab) = 1, equivalently the finite statistic 9-2 rank(D_tab)=7.

This is a truncated, untyped, distinct-support computation. It neither solves the
full Maurer--Cartan equation nor constructs a Floer differential. Every displayed
candidate is first checked for D_tab^2=0.
"""
import json
import itertools
from deform import rank_f2, square_entries_f2


def load(with_pent=True):
    with open("deform_full.json") as f:
        D = json.load(f)
    n = D["n"]
    d = D["d"]
    Tri = {int(a): D["Tri"][a] for a in D["Tri"]}
    Quad = {frozenset(int(x) for x in k.split(",")): v for k, v in D["Quad"].items()}
    Pent = {}
    if with_pent:
        try:
            with open("deform_pent.json") as f:
                P = json.load(f)
            Pent = {frozenset(int(x) for x in k.split(",")): v
                    for k, v in P["Pent"].items()}
        except FileNotFoundError:
            pass
    return n, d, Tri, Quad, Pent, D


def deformed(n, d, Tri, Quad, b, Pent=None):
    """Truncated table matrix for an untyped set of crossing indices.

    When Pent is supplied, include the table's distinct-support cubic term. This
    is not the all-order, ordered bounding-cochain deformation.
    """
    b = set(b)
    M = [row[:] for row in d]
    for S in b:
        T = Tri.get(S)
        if T:
            for i in range(n):
                for j in range(n):
                    M[i][j] ^= T[i][j]
    for key, Q in Quad.items():
        if key <= b:
            for i in range(n):
                for j in range(n):
                    M[i][j] ^= Q[i][j]
    if Pent:
        for key, Pm in Pent.items():
            if key <= b:
                for i in range(n):
                    for j in range(n):
                        M[i][j] ^= Pm[i][j]
    return M


def entries(n, M):
    return [(i, j) for i in range(n) for j in range(n) if M[i][j]]


if __name__ == "__main__":
    n, d, Tri, Quad, Pent, D = load()
    gP = D["gens"]
    order = "tri+quad+pent" if Pent else "tri+quad"
    raw_square = square_entries_f2(d)
    print(f"n={n}, bigons={D['bigons']}; finite polygon tables: {order} "
          f"({len(Pent)} pentagons)")
    print(f"D_big^2 audit: {'PASS' if not raw_square else 'FAIL'}; "
          f"nonzero entries {raw_square}")
    if raw_square:
        raise SystemExit("aborting before rank statistic: D_big^2 != 0")
    print(f"finite bigon statistic h_big=n-2 rank(D_big)="
          f"{n-2*rank_f2(d)} (rank(D_big)={rank_f2(d)})")

    act = set(S for S in Tri if any(any(r) for r in Tri[S]))
    for k in Quad:
        act |= set(k)
    for k in Pent:
        act |= set(k)
    active = sorted(act)
    print(f"active crossings: {len(active)} -> {active}")

    print("\n=== finite support screen ===")
    sols = []
    best = (99, None)
    maxr = min(len(active), 8)
    for r in range(1, maxr + 1):
        for combo in itertools.combinations(active, r):
            M = deformed(n, d, Tri, Quad, combo, Pent)
            rk = rank_f2(M)
            if rk < best[0]:
                best = (rk, combo)
            if rk == 1:
                sols.append(combo)
        if sols:
            break
        print(f"  size {r}: no support reaches the target statistic")
    if sols:
        audited = [(s, deformed(n, d, Tri, Quad, s, Pent)) for s in sols]
        failures = [(s, square_entries_f2(M)) for s, M in audited
                    if square_entries_f2(M)]
        square_zero = [(s, M) for s, M in audited if not square_entries_f2(M)]
        print(f"\nD_tab^2 audit of {len(sols)} candidate(s): "
              f"{len(square_zero)} PASS, {len(failures)} FAIL")
        for s, sq in failures:
            print(f"  FAIL support={s}: nonzero D_tab^2 entries {sq}")
        print(f"FINITE TARGET CANDIDATES ({len(sols)}), minimal size; "
              f"square-zero survivors {len(square_zero)}:")
        for s, M in square_zero[:40]:
            print(f"  support={s}; D_tab^2=0; entries {entries(n, M)}; "
                  f"h_tab={n-2*rank_f2(M)}")
    else:
        print(f"\nNO target-statistic subset up to size {maxr}.")
        if best[1]:
            best_matrix = deformed(n, d, Tri, Quad, best[1], Pent)
            best_square = square_entries_f2(best_matrix)
            print(f"  best-screen D_tab^2 audit: "
                  f"{'PASS' if not best_square else 'FAIL'}; "
                  f"nonzero entries {best_square}")
            if not best_square:
                print(f"  best square-zero support={best[1]}; "
                      f"h_tab={n-2*rank_f2(best_matrix)}; "
                      f"entries {entries(n, best_matrix)}")
