"""Twisted complexes over the algebra B of Kotelskiy, Watson and Zibrowius, and their morphism homology.

Conventions (Section 2.3 and Section 3.3 of the paper).

* An object is a pair (gens, arrows).  gens[i] is the idempotent of generator i: 'b' for iota_bullet and 'w' for
  iota_circ.  An arrow (i, kind, k, j) with kind in {'D', 'S'} and k >= 1 is the matrix entry kind^k from generator
  i to generator j.  In the paper's notation, D^k between iota_bullet generators is R^k, D^k between iota_circ
  generators is L^k, and S^k is M^k.  S^k joins generators of equal idempotent exactly when k is even.
* Products: same-face words concatenate, words in different faces multiply to zero, idempotents are units.
* A finite twisted complex is reduced if no entry of its differential is an idempotent.  All objects used in the
  paper are reduced.

Morphism homology.  For finite twisted complexes X and Y the morphism complex Mor(X, Y) is infinite dimensional.
Two exact methods are implemented, and they are independent of each other.

* hom_dim(X, Y): the element H = D + S^2 is central and B is a free F2[H]-module with basis 1, D on each diagonal
  block and S off the diagonal, so Mor(X, Y) is a finite free complex over F2[H].  From
  0 -> Mor -(H^n)-> Mor -> Mor/H^n -> 0 one gets dim H(Mor/H^n) = sum_i 2 min(n, n_i) + r n when
  H(Mor) = F2[H]^r + sum_i F2[H]/(H^{n_i}).  The free rank r is dim H(Mor at H = 1); the exponents n_i are read
  off from dim H(Mor/H^n), n = 1, 2, ..., which stops as soon as the increment equals r.
* hom_dim_cone(X, Y): the mapping-cone reduction of Lemma 3.5.  Mor is the cone of the chain map beta from the
  identity part (dimension e) to the positive part, which splits into three free complexes over F2[U_X],
  X = L, M, R, with U_L = L, U_R = R and U_M = M^2.  With h_m = dim H(C/U^m C) for the positive part C one has
  h_{m+1} - h_m = f + 2 #{i : deg p_i > m}, where f is the free rank of H(C).  So h_{m+1} = h_m implies f = 0 and
  that U^m annihilates H(C); then tau = h_m / 2, the map H(C) -> H(C/U^m C) is injective, and r = rank beta_* is
  computed in H(C/U^m C).  The dimension is e + tau - 2r.

Strict isomorphism of loop-type complexes (every generator meets at most two arrows) is decided by comparing the
words of the components up to rotation and reversal.
"""
import re
from collections import Counter

# ----------------------------------------------------------------------------------------------- objects

def mul(a, b):
    """Product of algebra elements a, b in path order (first a, then b); a = (kind, k) with kind in
    {'1', 'D', 'S'}.  Returns None for zero."""
    if a[0] == '1':
        return b
    if b[0] == '1':
        return a
    if a[0] != b[0]:
        return None
    return (a[0], a[1] + b[1])


def d_squared(o):
    """The nonzero entries of delta^2 (empty list iff delta^2 = 0)."""
    gens, arrows = o
    out = {}
    by_source = {}
    for (i, kind, k, j) in arrows:
        by_source.setdefault(i, []).append((kind, k, j))
    for (i, kind, k, j) in arrows:
        for (kind2, k2, l) in by_source.get(j, ()):
            p = mul((kind, k), (kind2, k2))
            if p is not None:
                key = (i, p, l)
                out[key] = out.get(key, 0) ^ 1
    return [key for key, v in out.items() if v]


def is_mc(o):
    """True iff delta^2 = 0, i.e. the arrows define a twisted complex (the Maurer--Cartan equation over B)."""
    return not d_squared(o)


def check_obj(o):
    gens, arrows = o
    for (i, kind, k, j) in arrows:
        assert kind in 'DS' and k >= 1, (i, kind, k, j)
        if kind == 'D':
            assert gens[i] == gens[j], ('D between different idempotents', i, j)
        else:
            assert (gens[i] == gens[j]) == (k % 2 == 0), ('parity of S', i, j, k)
    assert is_mc(o), 'delta^2 != 0'
    return o


def xor(arrows, *mods):
    """Symmetric difference of arrow lists (addition over F2)."""
    s = {}
    for t in list(arrows) + [x for m in mods for x in m]:
        s[t] = s.get(t, 0) ^ 1
    return [t for t, v in s.items() if v]


def deform(o, *terms):
    """The object (gens, delta + sum of terms) for elements of End(o) given as arrow lists."""
    return (list(o[0]), xor(o[1], *terms))


def union(*objs):
    gens, arrows = [], []
    for g, a in objs:
        off = len(gens)
        gens += list(g)
        arrows += [(i + off, kd, k, j + off) for (i, kd, k, j) in a]
    return (gens, arrows)


def components(o):
    """Connected components as objects (generators renumbered in increasing order)."""
    gens, arrows = o
    adj = {i: set() for i in range(len(gens))}
    for (i, _, _, j) in arrows:
        adj[i].add(j)
        adj[j].add(i)
    seen, comps = set(), []
    for v in range(len(gens)):
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        comp = []
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comp.sort()
        m = {x: t for t, x in enumerate(comp)}
        comps.append(([gens[x] for x in comp], [(m[i], kd, k, m[j]) for (i, kd, k, j) in arrows if i in m]))
    return comps


def valence(o):
    deg = Counter()
    for (i, _, _, j) in o[1]:
        deg[i] += 1
        deg[j] += 1
    return deg


def is_loop_type(o):
    """Every generator meets at most two arrows."""
    return all(v <= 2 for v in valence(o).values())


def arc_ends(o):
    deg = valence(o)
    return [v for v in range(len(o[0])) if deg[v] == 1]


# ----------------------------------------------------------------------------------------------- the paper's notation

def from_paper(labels, circ, arrows_text):
    """Build an object from the paper's notation.

    labels: generator labels in order (native index = position); circ: labels carrying iota_circ;
    arrows_text: 'a X b; ...' with X in {L, R, M} followed by an optional exponent, e.g. '13 M2 10'.
    Returns (object, labels)."""
    m = {v: t for t, v in enumerate(labels)}
    gens = ['w' if v in circ else 'b' for v in labels]
    return check_obj((gens, parse_arrows(arrows_text, m))), list(labels)


def parse_arrows(text, index):
    out = []
    for item in text.replace('\n', ' ').split(';'):
        item = item.strip()
        if not item:
            continue
        a, w, b = item.split()
        face = w[0]
        k = int(w[1:]) if len(w) > 1 else 1
        out.append((index[int(a)], 'S' if face == 'M' else 'D', k, index[int(b)]))
    return out


def label_of(o, t):
    """The paper's label of an arrow: R^k, L^k or M^k."""
    i, kind, k, j = t
    face = 'M' if kind == 'S' else ('R' if o[0][i] == 'b' else 'L')
    return face + (str(k) if k > 1 else '')


def arrow_str(o, t, labels=None):
    i, kind, k, j = t
    name = (lambda v: labels[v]) if labels is not None else (lambda v: v)
    return f'{name(i)} -{label_of(o, t)}-> {name(j)}'


def term_str(o, terms, labels=None):
    return ' + '.join(arrow_str(o, t, labels) for t in terms)


# ----------------------------------------------------------------------------------------------- text files

def read_complexes(path):
    """Read twisted complexes stored in the paper's notation:

        complex NAME
        gens 0 1 2 ...
        circ 3 4 ...
        arrows 1 M2 0; 2 R 1; ...
        end

    Returns {NAME: (object, labels)}."""
    out = {}
    cur = None
    for raw in open(path, encoding='utf-8'):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        key, _, rest = line.partition(' ')
        if key == 'complex':
            cur = {'name': rest.strip(), 'gens': [], 'circ': set(), 'arrows': ''}
        elif key == 'gens':
            cur['gens'] += [int(x) for x in rest.split()]
        elif key == 'circ':
            cur['circ'] |= {int(x) for x in rest.split()}
        elif key == 'arrows':
            cur['arrows'] += ' ' + rest + ';'
        elif key == 'end':
            out[cur['name']] = from_paper(cur['gens'], cur['circ'], cur['arrows'])
            cur = None
        else:
            raise ValueError('unexpected line: ' + line)
    return out


def write_complex(fh, name, o, labels=None):
    labels = labels if labels is not None else list(range(len(o[0])))
    fh.write(f'complex {name}\n')
    fh.write('gens ' + ' '.join(str(v) for v in labels) + '\n')
    fh.write('circ ' + ' '.join(str(labels[i]) for i in range(len(o[0])) if o[0][i] == 'w') + '\n')
    items = []
    for t in o[1]:
        i, kind, k, j = t
        face = 'M' if kind == 'S' else ('R' if o[0][i] == 'b' else 'L')
        items.append(f'{labels[i]} {face}{k if k > 1 else ""} {labels[j]}')
    for s in range(0, len(items), 8):
        fh.write('arrows ' + '; '.join(items[s:s + 8]) + '\n')
    fh.write('end\n')


# ----------------------------------------------------------------------------------------------- kht++ output

_GCH = {'⬮': 'b', '⬯': 'w'}     # black and white ellipses printed by kht++


def parse_kht_line(line):
    """Parse one component line of a kht++ output file (cxBNr-c2 or cxKhr-c2)."""
    s = line.strip()
    s = re.sub(r'^\d+\)\s*', '', s)
    s = re.sub(r'^h\^\s*-?\d+\s*q\^\s*-?[\d/]+\s*δ\^\s*-?[\d/]+\s*', '', s)
    s = re.split(r'\s{3,}', s)[0].strip()
    gens, seps, cur = [], [], ''
    for ch in s:
        if ch in _GCH:
            seps.append(cur)
            cur = ''
            gens.append(_GCH[ch])
        else:
            cur += ch
    seps.append(cur)
    assert seps[0] == '', ('text before the first generator', s)
    arrows = []
    closed = seps[-1].strip() != ''
    n = len(gens)
    for t in range(1, len(seps)):
        a = seps[t]
        if t == len(seps) - 1 and not closed:
            break
        i, j = t - 1, t % n
        fwd, bwd = '>' in a, '<' in a
        assert fwd != bwd, ('arrow direction', a)
        kind = 'S' if '~' in a else 'D'
        m = re.search(r'([DS])(?:\^(\d+))?', a)
        if m:
            assert m.group(1) == kind, a
            k = int(m.group(2)) if m.group(2) else 1
        else:
            k = 1 if kind == 'D' else (2 if gens[i] == gens[j] else 1)
        arrows.append((i, kind, k, j) if fwd else (j, kind, k, i))
    return check_obj((gens, arrows)), closed


def read_kht(path):
    """All components of a kht++ output file, as objects."""
    comps = []
    for line in open(path, encoding='utf-8'):
        if not line.strip() or line.lstrip().startswith('%'):
            continue
        comps.append(parse_kht_line(line)[0])
    return comps


# ----------------------------------------------------------------------------------------------- words and isomorphism

def walk(o):
    """Components of a loop-type object as (generator sequence, arrow labels, closed); labels[t] = (kind, k, dir)
    of the arrow between seq[t] and seq[t+1] (dir = +1 if it points forwards); for a closed component the last label
    joins seq[-1] to seq[0]."""
    g, a = o
    nb = {i: [] for i in range(len(g))}
    for t, (i, kd, k, j) in enumerate(a):
        nb[i].append((t, j, (kd, k, +1)))
        nb[j].append((t, i, (kd, k, -1)))
    assert all(len(v) <= 2 for v in nb.values()), 'not of loop type'
    seen, comps = set(), []
    order = [v for v in nb if len(nb[v]) <= 1] + list(range(len(g)))
    for start in order:
        if start in seen:
            continue
        seq, labs, used, cur, closed = [start], [], set(), start, False
        seen.add(start)
        while True:
            nxt = [e for e in nb[cur] if e[0] not in used]
            if not nxt:
                break
            t, y, lab = nxt[0]
            used.add(t)
            if y == start:
                labs.append(lab)
                closed = True
                break
            seq.append(y)
            labs.append(lab)
            seen.add(y)
            cur = y
        comps.append((seq, labs, closed))
    return comps


def _flip(lab):
    return (lab[0], lab[1], -lab[2])


def canon(o):
    """A complete invariant of strict isomorphism for loop-type objects: the sorted tuple of the least readings of
    the component words (arcs: two readings; closed components: every rotation and both orientations)."""
    g = o[0]
    out = []
    for seq, labs, closed in walk(o):
        if not closed:
            w1 = tuple([g[seq[0]]] + [x for t in range(len(labs)) for x in (labs[t], g[seq[t + 1]])])
            rs = seq[::-1]
            rl = [_flip(l) for l in reversed(labs)]
            w2 = tuple([g[rs[0]]] + [x for t in range(len(rl)) for x in (rl[t], g[rs[t + 1]])])
            out.append(('a', min(w1, w2)))
        else:
            n = len(seq)
            best = None
            for orient in (1, -1):
                if orient == 1:
                    s_, l_ = seq, labs
                else:
                    s_ = seq[::-1]
                    l_ = [_flip(labs[(n - 2 - t) % n]) for t in range(n)]
                for r in range(n):
                    w = tuple(x for t in range(n) for x in (g[s_[(r + t) % n]], l_[(r + t) % n]))
                    if best is None or w < best:
                        best = w
            out.append(('c', best))
    return tuple(sorted(out))


def strictly_isomorphic(o1, o2):
    """Strict isomorphism (an idempotent-preserving bijection of generators carrying one differential to the other
    entry by entry) of loop-type objects."""
    if len(o1[0]) != len(o2[0]) or len(o1[1]) != len(o2[1]):
        return False
    if not (is_loop_type(o1) and is_loop_type(o2)):
        raise ValueError('strict isomorphism is implemented for loop-type objects only')
    return canon(o1) == canon(o2)


def component_sizes(o):
    return sorted(len(c[0]) for c in components(o))


def word_str(o):
    """Human-readable word of a connected loop-type object (o = iota_circ, * = iota_bullet)."""
    (seq, labs, closed), = walk(o)
    g = o[0]
    s = 'o' if g[seq[0]] == 'w' else '*'
    for t, lab in enumerate(labs):
        kd, k, dr = lab
        name = kd + (str(k) if k > 1 else '')
        s += ('-%s->' % name) if dr > 0 else ('<-%s-' % name)
        y = seq[(t + 1) % len(seq)]
        s += 'o' if g[y] == 'w' else '*'
    return s + (' (closed)' if closed else '')


# ----------------------------------------------------------------------------------------------- F2 linear algebra

def _rank(vectors):
    pivots = {}
    r = 0
    for v in vectors:
        while v:
            low = v & -v
            p = pivots.get(low)
            if p is None:
                pivots[low] = v
                r += 1
                break
            v ^= p
    return r


# ----------------------------------------------------------------------------------------------- method 1: F2[H]

def _to_H(gx, gy, elt):
    """Coordinates of the algebra element elt from idempotent gx to gy in the F2[H]-basis {1, D} (same idempotent)
    or {S} (opposite idempotents): a list of (basis symbol, power of H)."""
    kind, k = elt
    if kind == '1':
        assert gx == gy
        return [('1', 0)]
    if kind == 'D':
        assert gx == gy
        return [('D', k - 1)]
    if gx == gy:                       # S^{2m} = H^m 1 + H^{m-1} D
        assert k % 2 == 0
        m = k // 2
        return [('1', m), ('D', m - 1)]
    assert k % 2 == 1
    return [('S', (k - 1) // 2)]


_BASIS_ELT = {'1': ('1', 0), 'D': ('D', 1), 'S': ('S', 1)}


def mor_complex_H(X, Y):
    """Mor(X, Y) as a free F2[H]-complex: (basis, columns), columns[c] = list of (row, power of H)."""
    gx, ax = X
    gy, ay = Y
    basis = []
    for x in range(len(gx)):
        for y in range(len(gy)):
            if gx[x] == gy[y]:
                basis.append((x, '1', y))
                basis.append((x, 'D', y))
            else:
                basis.append((x, 'S', y))
    idx = {b: t for t, b in enumerate(basis)}
    into_x = {}
    for (i, kd, k, j) in ax:
        into_x.setdefault(j, []).append(((kd, k), i))
    out_y = {}
    for (i, kd, k, j) in ay:
        out_y.setdefault(i, []).append(((kd, k), j))
    cols = []
    for (x, sym, y) in basis:
        a = _BASIS_ELT[sym]
        acc = {}
        for (b, z) in out_y.get(y, ()):
            p = mul(a, b)
            if p is None:
                continue
            for (s2, pw) in _to_H(gx[x], gy[z], p):
                key = (idx[(x, s2, z)], pw)
                acc[key] = acc.get(key, 0) ^ 1
        for (c, x2) in into_x.get(x, ()):
            p = mul(c, a)
            if p is None:
                continue
            for (s2, pw) in _to_H(gx[x2], gy[y], p):
                key = (idx[(x2, s2, y)], pw)
                acc[key] = acc.get(key, 0) ^ 1
        cols.append([kk for kk, v in acc.items() if v])
    return basis, cols


def _quotient_dim(basis, cols, n):
    vecs = []
    for col in cols:
        for p in range(n):
            v = 0
            for (row, pw) in col:
                q = p + pw
                if q < n:
                    v ^= 1 << (row * n + q)
            vecs.append(v)
    return len(basis) * n - 2 * _rank(vecs)


def _free_rank(basis, cols):
    vecs = []
    for col in cols:
        v = 0
        for (row, pw) in col:
            v ^= 1 << row
        vecs.append(v)
    return len(basis) - 2 * _rank(vecs)


def mor_module(X, Y, nmax=60):
    """(r, sorted exponents n_i, [d(1), d(2), ...]) with H(Mor(X, Y)) = F2[H]^r + sum_i F2[H]/(H^{n_i})."""
    basis, cols = mor_complex_H(X, Y)
    r = _free_rank(basis, cols)
    d = [0]
    n = 0
    while True:
        n += 1
        d.append(_quotient_dim(basis, cols, n))
        inc = d[n] - d[n - 1]
        assert inc >= r and (inc - r) % 2 == 0, (inc, r, d)
        if inc == r:
            break
        if n >= nmax:
            raise RuntimeError('torsion exponent > %d' % nmax)
    incs = [d[t] - d[t - 1] for t in range(1, len(d))]
    ns = []
    for t in range(1, len(incs) + 1):
        cnt = (incs[t - 1] - r) // 2
        nxt = (incs[t] - r) // 2 if t < len(incs) else 0
        ns += [t] * (cnt - nxt)
    return r, sorted(ns), d[1:]


def hom_dim(X, Y):
    """dim_F2 H Mor(X, Y), by the F2[H]-module structure; raises if the homology is infinite."""
    r, ns, d = mor_module(X, Y)
    if r:
        raise ValueError('H Mor(X, Y) has a free F2[H]-summand of rank %d (infinite dimensional)' % r)
    return sum(ns)


# ----------------------------------------------------------------------------------------------- method 2: Lemma 3.5

def _face(kind, idem):
    if kind == 'S':
        return 'M'
    return 'R' if idem == 'b' else 'L'


def hom_dim_cone(X, Y, mmax=40):
    """dim H Mor(X, Y) by the mapping-cone reduction of Lemma 3.5.  Returns a dict with e, tau, r, f (free ranks of
    the three face complexes), m (a power of U that annihilates H of the positive part) and dim.  dim is None if some
    f_X > 0 (infinite homology)."""
    gx, ax = X
    gy, ay = Y
    into_x = {}
    for (i, kd, k, j) in ax:
        into_x.setdefault(j, []).append(((kd, k), i))
    out_y = {}
    for (i, kd, k, j) in ay:
        out_y.setdefault(i, []).append(((kd, k), j))
    # positive part: free generators (x, y, face) with base exponent k0 and step
    gens_pos = []
    for x in range(len(gx)):
        for y in range(len(gy)):
            if gx[x] == gy[y]:
                gens_pos.append((x, y, 'D', 1, 1))                     # face L or R
                gens_pos.append((x, y, 'S', 2, 2))                     # face M, even powers
            else:
                gens_pos.append((x, y, 'S', 1, 2))                     # face M, odd powers
    gidx = {(x, y, kd): t for t, (x, y, kd, k0, st) in enumerate(gens_pos)}

    def d_term(x, a, y):
        """d of the entry (x, a, y) as a list of (x', (kind, k), y')."""
        out = []
        for (b, z) in out_y.get(y, ()):
            p = mul(a, b)
            if p is not None:
                out.append((x, p, z))
        for (c, x2) in into_x.get(x, ()):
            p = mul(c, a)
            if p is not None:
                out.append((x2, p, y))
        return out

    def coords(term):
        """(generator index, U-power) of a positive entry."""
        x, (kd, k), y = term
        t = gidx[(x, y, kd)]
        k0, st = gens_pos[t][3], gens_pos[t][4]
        assert k >= k0 and (k - k0) % st == 0
        return t, (k - k0) // st

    # columns over F2[U]: for each free generator, d(generator) as (row, U-power) pairs
    cols = []
    for (x, y, kd, k0, st) in gens_pos:
        acc = {}
        for term in d_term(x, (kd, k0), y):
            key = coords(term)
            acc[key] = acc.get(key, 0) ^ 1
        cols.append([kk for kk, v in acc.items() if v])
    faces = [_face(kd, gx[x]) for (x, y, kd, k0, st) in gens_pos]
    # free ranks of the face complexes: H at U = 1
    f = {}
    for F in 'LMR':
        rows = [t for t in range(len(gens_pos)) if faces[t] == F]
        pos = {t: s for s, t in enumerate(rows)}
        vecs = []
        for t in rows:
            v = 0
            for (row, pw) in cols[t]:
                assert faces[row] == F
                v ^= 1 << pos[row]
            vecs.append(v)
        f[F] = len(rows) - 2 * _rank(vecs)
    # identity part and beta
    ident = [(x, y) for x in range(len(gx)) for y in range(len(gy)) if gx[x] == gy[y]]
    e = len(ident)
    beta = []
    for (x, y) in ident:
        acc = {}
        for term in d_term(x, ('1', 0), y):
            key = coords(term)
            acc[key] = acc.get(key, 0) ^ 1
        beta.append([kk for kk, v in acc.items() if v])
    res = {'e': e, 'f': f}
    if any(f.values()):
        res.update(tau=None, r=None, m=None, dim=None)
        return res
    n = len(gens_pos)

    def trunc_vecs(m):
        vecs = []
        for col in cols:
            for p in range(m):
                v = 0
                for (row, pw) in col:
                    q = p + pw
                    if q < m:
                        v ^= 1 << (row * m + q)
                vecs.append(v)
        return vecs

    h_prev, m = None, 0
    while True:
        m += 1
        vecs = trunc_vecs(m)
        rk = _rank(vecs)
        h = n * m - 2 * rk
        if h_prev is not None and h == h_prev:
            m -= 1                       # U^m annihilates H(C) for the previous m
            break
        h_prev = h
        if m >= mmax:
            raise RuntimeError('torsion exponent > %d' % mmax)
    vecs = trunc_vecs(m)
    rk = _rank(vecs)
    tau = (n * m - 2 * rk) // 2
    bvecs = []
    for col in beta:
        v = 0
        for (row, pw) in col:
            if pw < m:
                v ^= 1 << (row * m + pw)
        bvecs.append(v)
    r = _rank(vecs + bvecs) - rk
    res.update(tau=tau, r=r, m=m, dim=e + tau - 2 * r)
    return res


# ----------------------------------------------------------------------------------------------- reporting

class Report:
    """Prints 'PASS'/'FAIL' lines with the paper's value next to the computed one, and keeps count."""

    def __init__(self, title):
        self.fail = 0
        self.npass = 0
        print(title)
        print('=' * len(title))

    def check(self, what, paper, computed, source='paper'):
        """source = 'paper' for a value stated in the paper, 'expected' for an additional consistency check."""
        ok = paper == computed
        self.npass += ok
        self.fail += not ok
        print(f'{"PASS" if ok else "FAIL"}  {what}: {source} {paper}, computed {computed}', flush=True)
        return ok

    def note(self, text):
        print('      ' + text, flush=True)

    def section(self, text):
        print('\n-- ' + text, flush=True)

    def finish(self, t0=None):
        import time
        extra = f' ({time.time() - t0:.0f} s)' if t0 is not None else ''
        print(f'\n{self.npass} checks passed, {self.fail} failed{extra}')
        return 1 if self.fail else 0
