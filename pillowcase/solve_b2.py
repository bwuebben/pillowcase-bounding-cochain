#!/usr/bin/env python3
"""
solve_b2.py -- step (e), part 5: search for b_2 (RESEARCH_LOG sec 29/30).

Loads the cached deformed-differential tables (deform_full.json: d, Tri per
crossing, Quad per crossing-pair) and evaluates, for any candidate support
b subset {crossings},

    partial_b[i][j] = d[i][j]  XOR  (+)_{S in b} Tri[S][i][j]
                               XOR  (+)_{ {S,S'} subset b } Quad[{S,S'}][i][j]

then searches all subsets of the ACTIVE crossings (those appearing in Tri or Quad)
for rank(partial_b) = 1  <=>  rank HF = 9 - 2 = 7 = rank I^natural(P(-2,3,5)).

This is the ground-truth check on the sec-29 obstruction: if no tri+quad subset
gives rank 1, the pentagon (mu^4) terms are genuinely required.
"""
import json
import itertools
from deform import rank_f2


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
    """partial_b for support b (an iterable of crossing indices), including the
    cubic mu^4 pentagon term when Pent is supplied."""
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
    print(f"n={n}, bigons={D['bigons']}, undeformed rank={rank_f2(d)} "
          f"HF={n-2*rank_f2(d)}; polygon order: {order} ({len(Pent)} pentagons)")

    act = set(S for S in Tri if any(any(r) for r in Tri[S]))
    for k in Quad:
        act |= set(k)
    for k in Pent:
        act |= set(k)
    active = sorted(act)
    print(f"active crossings: {len(active)} -> {active}")

    print(f"\n=== searching subsets of active crossings for rank(partial_b)=1 ===")
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
            print(f"  size {r}: found {len(sols)} rank-1 support(s)")
            break
        print(f"  size {r}: none (running best rank {best[0]} at {best[1]})")
    if sols:
        print(f"\nRANK-1 SOLUTIONS ({len(sols)}), minimal size:")
        for s in sols[:40]:
            M = deformed(n, d, Tri, Quad, s, Pent)
            print(f"  b = {s}  partial_b entries {entries(n, M)}")
    else:
        print(f"\nNO rank-1 subset up to size {maxr}. Best rank {best[0]} at b={best[1]}")
        if best[1]:
            print(f"  entries: {entries(n, deformed(n, d, Tri, Quad, best[1], Pent))}")
