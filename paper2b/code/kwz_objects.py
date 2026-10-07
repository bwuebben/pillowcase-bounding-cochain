"""The explicit twisted complexes of Sections 3 and 4 and Appendix A, in the paper's notation, and loaders for the
data files.

Generators carry the labels printed in the paper.  N_7 is the complex N of Section 3.3, E and E_* are the earring
complexes (3.2) and (3.3) of the first and second closures, N_5 and b_5 are (A.1) and (A.2), and the corner terms
kappa_q and kappa'_q are those of Theorem 4.1(c) and Remark 4.3.
"""
import os

import kwz_algebra as A

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

# ---------------------------------------------------------------- Section 3.3
N7_LABELS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 26, 27, 30, 31, 32, 33,
             34, 35, 36]
N7_CIRC = {1, 2, 5, 6, 15, 16, 20, 21, 32, 33}
N7_ARROWS = """17 R 18; 16 M 17; 16 L 15; 15 M 14; 14 R 13; 13 M2 10; 10 R 9; 8 M2 9; 7 R 8; 6 M 7; 6 L 5;
5 M 4; 4 R 3; 2 M 3; 2 L 1; 1 M 0; 0 R 19; 20 M 19; 20 L 21; 21 M 22; 22 R 23; 23 M2 26; 26 R 27;
30 M2 27; 31 R 30; 32 M 31; 32 L 33; 33 M 34; 34 R 35; 35 M2 36"""

E_LABELS = [0, 3, 4, 5, 6, 7]
E_CIRC = {5, 6}
E_ARROWS = "3 M2 0; 4 R 3; 5 M 4; 6 L 5; 6 M 7; 0 R 7"

ESTAR_LABELS = list(range(14))
ESTAR_CIRC = {1, 2, 7, 8, 11, 12}
ESTAR_ARROWS = """1 M 0; 1 L 2; 2 M 3; 3 R 4; 5 M2 4; 6 R 5; 7 M 6; 8 L 7; 8 M 9; 10 R 9; 11 M 10; 12 L 11;
12 M 13; 0 R 13"""

# the four switches of the table in Section 3.3: (old arrows, new arrows)
SWITCHES7 = {'S18': ("4 R 3; 22 R 23", "22 R 3; 4 R 23"),
             'S25': ("6 L 5; 20 L 21", "20 L 5; 6 L 21"),
             'S69': ("0 R 19; 26 R 27", "0 R 27; 26 R 19"),
             'S74': ("20 L 21; 32 L 33", "20 L 33; 32 L 21")}

# ---------------------------------------------------------------- Appendix A.1
N5_LABELS = list(range(23))
N5_CIRC = {3, 4, 9, 10, 13, 14, 19, 20}
N5_ARROWS = """1 M2 0; 2 R 1; 3 M 2; 4 L 3; 4 M 5; 5 R 6; 6 M2 7; 8 R 7; 9 M 8; 10 L 9; 10 M 11; 12 R 11;
13 M 12; 14 L 13; 14 M 15; 16 R 15; 17 M2 16; 18 R 17; 19 M 18; 20 L 19; 20 M 21; 21 R 22"""
B5 = "8 R 7; 16 R 15; 8 R 15; 16 R 7"

# ---------------------------------------------------------------- corner terms (Theorem 4.1(c), Remark 4.3)
KAPPA = {7: "0 R 19; 36 R 19", 5: "12 R 11; 0 R 11"}
KAPPA_PRIME = {7: "2 M 3; 18 M2 3", 5: "10 M 11; 22 M2 11"}


def _terms(text, labels):
    idx = {v: t for t, v in enumerate(labels)}
    return A.parse_arrows(text, idx)


def N7():
    return A.from_paper(N7_LABELS, N7_CIRC, N7_ARROWS)


def E():
    return A.from_paper(E_LABELS, E_CIRC, E_ARROWS)


def Estar():
    return A.from_paper(ESTAR_LABELS, ESTAR_CIRC, ESTAR_ARROWS)


def N5():
    return A.from_paper(N5_LABELS, N5_CIRC, N5_ARROWS)


def N(q):
    """(N_q, labels) for q = 5, 7."""
    return {5: N5, 7: N7}[q]()


def switch7(name):
    """b_s in End(N_7), as an arrow list (old arrows + new arrows)."""
    old, new = SWITCHES7[name]
    return _terms(old + ';' + new, N7_LABELS)


def b5():
    return _terms(B5, N5_LABELS)


def selected_b(q):
    """The selected Maurer--Cartan elements: {name: arrow list} (b_{S_18}, b_{S_25} at q = 7; b_5 at q = 5)."""
    if q == 7:
        return {'b_S18': switch7('S18'), 'b_S25': switch7('S25')}
    return {'b_5': b5()}


def kappa(q):
    return _terms(KAPPA[q], N(q)[1])


def kappa_prime(q):
    return _terms(KAPPA_PRIME[q], N(q)[1])


def end_labels(q):
    """The labels of the ends of N_q at (0, 0) and at (pi, 0)."""
    return {7: {'(0,0)': 36, '(pi,0)': 18}, 5: {'(0,0)': 0, '(pi,0)': 22}}[q]


# ---------------------------------------------------------------- data files

def complexes(name):
    """Twisted complexes stored in data/<name>.txt: {name: (object, labels)}."""
    return A.read_complexes(os.path.join(DATA, name + '.txt'))


def kht_bn(q):
    """KWZ's BN~(T_q) as computed by kht++ (data/kht/T<q>yy/cxBNr-c2)."""
    comps = A.read_kht(os.path.join(DATA, 'kht', f'T{q}yy', 'cxBNr-c2'))
    assert len(comps) == 1
    return comps[0]


def kht_kh(q):
    """The components of KWZ's Kh~(T_q) as computed by kht++."""
    return A.read_kht(os.path.join(DATA, 'kht', f'T{q}yy', 'cxKhr-c2'))


def khovanov_ranks():
    """{(slope string, q): rank of reduced Khovanov homology of K_r(q) = num(Q_r + Q_{1/3} + Q_{1/q}) over F2},
    from data/khovanov_ranks.txt; slope 'inf' is r = infinity."""
    out = {}
    for line in open(os.path.join(DATA, 'khovanov_ranks.txt'), encoding='utf-8'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        q, r, rk = line.split()
        out[(r, int(q))] = int(rk)
    return out


def slope_key(r):
    return 'inf' if r is None else str(r)


def slope_set(q):
    """The slopes r of Corollary 4.6 and the paragraph after it: every slope with |r| <= 2 and denominator at most
    14, every slope with |r| <= 5 and denominator at most 6, the integers with |r| <= 10, +-r_q, and infinity (None)."""
    import math
    from fractions import Fraction as F
    S = set()
    for d in range(1, 15):
        for p in range(-2 * d, 2 * d + 1):
            if math.gcd(p, d) == 1:
                S.add(F(p, d))
    for d in range(1, 7):
        for p in range(-5 * d, 5 * d + 1):
            if math.gcd(p, d) == 1:
                S.add(F(p, d))
    for p in range(6, 11):
        S.add(F(p))
        S.add(F(-p))
    S.add(F(q - 1, 3 * q - 4))
    S.add(F(-(q - 1), 3 * q - 4))
    return sorted(S) + [None]


def r_q(q):
    from fractions import Fraction as F
    return F(1 - q, 3 * q - 4)


def read_kht_rational():
    """{slope string: (BN~(Q_s), Kh~(Q_s))} from data/kht_rational.txt (kht++ output; the slope is kht++'s)."""
    out = {}
    cur = None
    for raw in open(os.path.join(DATA, 'kht_rational.txt'), encoding='utf-8'):
        line = raw.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        key, _, rest = line.partition(' ')
        if key == 'slope':
            cur = {'s': rest.strip(), 'bn': [], 'kh': None}
        elif key == 'word':
            cur['word'] = rest.strip()
        elif key == 'BN':
            cur['bn'].append(A.parse_kht_line(rest)[0])
        elif key == 'Kh':
            cur['kh'] = A.parse_kht_line(rest)[0]
        elif key == 'end':
            out[cur['s']] = (A.union(*cur['bn']), cur['kh'])
    return out
