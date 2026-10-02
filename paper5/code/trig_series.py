"""Exact arithmetic for trigonometric polynomials in xi whose frequencies are linear in h.

A term is  c(h) * cos(w(h) xi)  or  c(h) * sin(w(h) xi),  with c a polynomial in h with rational
coefficients and w(h) = alpha + beta h, alpha and beta rational.  Products are expanded by the
product-to-sum formulas, so every expression built from such terms is again a finite sum of them.
This is all that is needed to compute Taylor coefficients in xi (as polynomials in h) exactly, and to
bound derivatives in xi by sums |c| |w|^k in interval arithmetic.
"""
from fractions import Fraction as Fr
from math import factorial


# ----------------------------------------------------------------------------- polynomials in h
def padd(p, q):
    n = max(len(p), len(q))
    out = [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]
    return ptrim(out)


def pmul(p, q):
    if not p or not q:
        return []
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return ptrim(out)


def pscale(p, c):
    return ptrim([c * a for a in p])


def ptrim(p):
    p = [Fr(a) for a in p]
    while p and p[-1] == 0:
        p.pop()
    return p


def ppow(p, k):
    out = [Fr(1)]
    for _ in range(k):
        out = pmul(out, p)
    return out


def pderiv(p):
    return ptrim([i * p[i] for i in range(1, len(p))])


def pdiv_h(p):
    """p(h)/h; p must vanish at h = 0."""
    p = ptrim(p)
    if p and p[0] != 0:
        raise ValueError("polynomial does not vanish at h = 0")
    return p[1:]


def pfrom_roots(lead, roots):
    """lead * prod (h - r)."""
    out = [Fr(lead)]
    for r in roots:
        out = pmul(out, [Fr(-r), Fr(1)])
    return out


def pstr(p):
    out = ""
    for i, a in enumerate(p):
        if not a:
            continue
        mono = "" if i == 0 else ("*h" if i == 1 else f"*h^{i}")
        sgn = "-" if a < 0 else "+"
        out += (f" {sgn} " if out else ("-" if a < 0 else "")) + f"{abs(a)}{mono}"
    return out or "0"


# ----------------------------------------------------------------------------- trigonometric polynomials
def _canon(kind, a, b):
    """normalize the frequency so that (a, b) > (0, 0) lexicographically; returns (kind, a, b, sign)."""
    if (a, b) < (0, 0):
        return kind, -a, -b, (1 if kind == 'c' else -1)
    return kind, a, b, 1


class TrigPoly:
    """sum over keys (kind, alpha, beta) of coef(h) * trig((alpha + beta h) xi), kind in {'c', 's'}."""

    def __init__(self, terms=None):
        self.t = {}
        for key, c in (terms or {}).items():
            self._add(key, c)

    def _add(self, key, c):
        kind, a, b = key
        a, b = Fr(a), Fr(b)
        if a == 0 and b == 0:
            if kind == 's':
                return
        kind, a, b, sg = _canon(kind, a, b)
        c = pscale(c, sg)
        key = (kind, a, b)
        new = padd(self.t.get(key, []), c)
        if new:
            self.t[key] = new
        elif key in self.t:
            del self.t[key]

    @staticmethod
    def const(c):
        return TrigPoly({('c', 0, 0): [Fr(c)] if not isinstance(c, list) else c})

    @staticmethod
    def cos(a, b, coef=(1,)):
        return TrigPoly({('c', a, b): [Fr(x) for x in coef]})

    @staticmethod
    def sin(a, b, coef=(1,)):
        return TrigPoly({('s', a, b): [Fr(x) for x in coef]})

    def __add__(self, o):
        out = TrigPoly(self.t)
        for k, c in o.t.items():
            out._add(k, c)
        return out

    def __neg__(self):
        return TrigPoly({k: pscale(c, -1) for k, c in self.t.items()})

    def __sub__(self, o):
        return self + (-o)

    def scale(self, p):
        """multiply by a polynomial in h (list) or a number."""
        p = p if isinstance(p, list) else [Fr(p)]
        return TrigPoly({k: pmul(c, p) for k, c in self.t.items()})

    def __mul__(self, o):
        if not isinstance(o, TrigPoly):
            return self.scale(o)
        out = TrigPoly()
        half = Fr(1, 2)
        for (k1, a1, b1), c1 in self.t.items():
            for (k2, a2, b2), c2 in o.t.items():
                c = pscale(pmul(c1, c2), half)
                s, d = (a1 + a2, b1 + b2), (a1 - a2, b1 - b2)
                if k1 == 'c' and k2 == 'c':
                    out._add(('c',) + d, c); out._add(('c',) + s, c)
                elif k1 == 's' and k2 == 's':
                    out._add(('c',) + d, c); out._add(('c',) + s, pscale(c, -1))
                elif k1 == 's' and k2 == 'c':
                    out._add(('s',) + s, c); out._add(('s',) + d, c)
                else:  # cos A sin B = (sin(A+B) - sin(A-B)) / 2
                    out._add(('s',) + s, c); out._add(('s',) + d, pscale(c, -1))
        return out

    __rmul__ = __mul__

    def at_h0(self):
        """the trigonometric polynomial in xi obtained by setting h = 0 (dict (kind, alpha) -> value)."""
        out = {}
        for (k, a, b), c in self.t.items():
            v = c[0] if c else Fr(0)
            key = (k, a)
            out[key] = out.get(key, Fr(0)) + v
        return {k: v for k, v in out.items() if v != 0 and not (k[0] == 's' and k[1] == 0)}

    def taylor(self, order):
        """Taylor coefficients in xi up to xi^order, each a polynomial in h."""
        coeffs = [[] for _ in range(order + 1)]
        for (k, a, b), c in self.t.items():
            w = [a, b]
            for j in range(order + 1):
                if k == 'c' and j % 2 == 1:
                    continue
                if k == 's' and j % 2 == 0:
                    continue
                sgn = (-1) ** (j // 2)
                coeffs[j] = padd(coeffs[j], pscale(pmul(c, ppow(w, j)), Fr(sgn, factorial(j))))
        return coeffs

    def dxi(self):
        """derivative in xi."""
        out = TrigPoly()
        for (k, a, b), c in self.t.items():
            cw = pmul(c, [a, b])
            if k == 'c':
                out._add(('s', a, b), pscale(cw, -1))
            else:
                out._add(('c', a, b), cw)
        return out

    def dh(self):
        """derivative in h, returned as a pair (T0, T1) with d/dh = T0 + xi * T1."""
        T0, T1 = TrigPoly(), TrigPoly()
        for (k, a, b), c in self.t.items():
            T0._add((k, a, b), pderiv(c))
            if b:
                cb = pscale(c, b)
                if k == 'c':
                    T1._add(('s', a, b), pscale(cb, -1))
                else:
                    T1._add(('c', a, b), cb)
        return T0, T1

    def items(self):
        return self.t.items()
