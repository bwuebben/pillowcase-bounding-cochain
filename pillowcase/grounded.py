#!/usr/bin/env python3
"""Grounded tangle character varieties from REP THEORY -- no Smith conventions.
A 2-string tangle's pillowcase curve = {traceless SU(2) reps of the tangle group,
restricted to the 4-punctured boundary} in (gamma,theta) coords. We compute the
boundary meridians X1..X4 as quaternion WORDS in the arc-meridian generators,
solve for traceless reps, and read off (gamma,theta). Ground truth, no Psi."""
import math
Q1=(1,0,0,0)
def qm(p,q):
    w1,x1,y1,z1=p; w2,x2,y2,z2=q
    return (w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
            w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2)
def qinv(q): w,x,y,z=q; return (w,-x,-y,-z)   # unit quaternion inverse
def tracefree(axis):  # pure-imaginary unit quaternion from a unit 3-vector
    x,y,z=axis; return (0.0,x,y,z)
def word(gens, w):    # w = list of (gen_index, +/-1)
    r=Q1
    for gi,s in w:
        g=gens[gi]; r=qm(r, g if s>0 else qinv(g))
    return r

def pillowcase_coords(Xs):
    """Given X1..X4 (traceless, product=1), return (gamma,theta) in the pillowcase,
    or None if reducible/degenerate. Uses: A=X1X2 has axis n; the Xi are coplanar
    perp to n; measure their angles beta_m in that plane; gamma=beta2-beta1, theta=beta4."""
    X1,X2,X3,X4=Xs
    A=qm(X1,X2)
    n=(A[1],A[2],A[3]); nn=math.sqrt(sum(c*c for c in n))
    if nn<1e-9: return None  # A central: reducible locus / corner
    n=tuple(c/nn for c in n)
    # build an orthonormal frame (e1,e2) in the plane perp to n, e1 = X1 direction
    v1=(X1[1],X1[2],X1[3])
    # project out n
    dot=sum(v1[i]*n[i] for i in range(3))
    e1=tuple(v1[i]-dot*n[i] for i in range(3)); e1n=math.sqrt(sum(c*c for c in e1))
    if e1n<1e-9: return None
    e1=tuple(c/e1n for c in e1)
    e2=(n[1]*e1[2]-n[2]*e1[1], n[2]*e1[0]-n[0]*e1[2], n[0]*e1[1]-n[1]*e1[0])  # n x e1
    def beta(X):
        v=(X[1],X[2],X[3])
        return math.atan2(sum(v[i]*e2[i] for i in range(3)), sum(v[i]*e1[i] for i in range(3)))
    b=[beta(X) for X in Xs]
    gamma=(b[1]-b[0])%(2*math.pi); theta=(b[3]-b[0])%(2*math.pi)
    return (gamma, theta)

def tangle_curve(boundary_words, ngen=2, N=240):
    """Trace the traceless character variety of a 2-string tangle whose boundary
    meridians X1..X4 are the given WORDS in `ngen` arc-meridian generators.
    Parametrize reps: generator 0 = i (fixed); generator 1 = axis at angle eta in
    the i-j plane; (for ngen>2, further generators need more params -- here ngen=2)."""
    pts=[]
    for k in range(N+1):
        eta=math.pi*k/N
        gens=[tracefree((1,0,0)),                      # m1 = i
              tracefree((math.cos(eta),math.sin(eta),0))]  # m2 at angle eta
        Xs=[word(gens,w) for w in boundary_words]
        # check boundary relation X1X2X3X4 = 1
        prod=Q1
        for X in Xs: prod=qm(prod,X)
        if abs(prod[0]-1)>1e-6: continue   # not a valid boundary rep
        c=pillowcase_coords(Xs)
        if c: pts.append(c)
    return pts

if __name__=="__main__":
    # TRIVIAL infinity-tangle: strands connect punctures 1-2 and 3-4, untwisted.
    #   X1=m1, X2=m1^{-1}? Try the version that satisfies X1X2X3X4=1 with tracefree gens.
    #   Two arcs, meridians m1 (arc through punctures 1,2) and m2 (through 3,4).
    #   Standard: X1=m1, X2=m1, X3=m2, X4=m2 gives X1X2X3X4=m1^2 m2^2 = (-1)(-1)=1. OK.
    inf_tangle=[[(0,1)],[(0,1)],[(1,1)],[(1,1)]]   # X1=X2=m1, X3=X4=m2
    zero_tangle=[[(0,1)],[(1,1)],[(1,1)],[(0,1)]]  # X1=m1,X2=m2,X3=m2,X4=m1 -> m1 m2 m2 m1=1
    for name,tw in [("infinity",inf_tangle),("zero",zero_tangle)]:
        pts=tangle_curve(tw)
        if not pts: print(f"{name}-tangle: (no valid reps)"); continue
        gam=[round(math.degrees(p[0])) for p in pts]; th=[round(math.degrees(p[1])) for p in pts]
        print(f"{name}-tangle: {len(pts)} pts; gamma in [{min(gam)},{max(gam)}] deg, "
              f"theta in [{min(th)},{max(th)}] deg  "
              f"=> {'gamma=0 (vertical edge)' if max(gam)<2 else ('theta=0 (horiz edge)' if max(th)<2 else 'diagonal')}")


# --- braid/twist action on the 4 boundary meridian WORDS (rational tangles) ---
def winv(w): return [(gi,-s) for gi,s in reversed(w)]
def braid(words, i, inv=False):
    """Artin generator sigma_i on boundary words (i in {0,1,2}). Builds rational
    tangles by twisting adjacent strands: sigma_i: W_i -> W_i W_{i+1} W_i^{-1},
    W_{i+1} -> W_i."""
    W=list(words)
    if not inv:
        Wi=W[i]; W[i]=Wi+W[i+1]+winv(Wi); W[i+1]=Wi
    else:
        Wi1=W[i+1]; W[i+1]=winv(Wi1)+W[i]+Wi1; W[i]=Wi1
    return W

def rational_tangle(twists):
    """Build a rational tangle from the 0-tangle by a sequence of twists.
    twists = list of (i, inv) braid generators."""
    W=[[(0,1)],[(1,1)],[(1,1)],[(0,1)]]   # 0-tangle
    for (i,inv) in twists: W=braid(W,i,inv)
    return W

if __name__=="__main__":
    print("\n--- rational (twisted) tangles: do they trace straight lines? ---")
    tests = {
      "1 twist (sigma_2)":        [(1,False)],
      "2 twists (sigma_2^2)":     [(1,False),(1,False)],
      "3 twists (sigma_2^3)":     [(1,False),(1,False),(1,False)],
      "mixed s2 s1":              [(1,False),(0,False)],
    }
    for name,tw in tests.items():
        W=rational_tangle(tw)
        pts=tangle_curve(W)
        if len(pts)<4: print(f"  {name:22s}: {len(pts)} pts (degenerate/edge)"); continue
        # fit a line theta = m*gamma + c on the interior points; report slope & straightness
        import statistics
        gs=[p[0] for p in pts]; ts=[p[1] for p in pts]
        # unwrap theta for a clean fit near the middle
        gm=[math.degrees(g) for g in gs]; tm=[math.degrees(t) for t in ts]
        # crude slope from endpoints of the sorted-by-gamma interior
        order=sorted(range(len(pts)), key=lambda k:gs[k])
        g0,g1=gs[order[len(order)//4]],gs[order[3*len(order)//4]]
        t0,t1=ts[order[len(order)//4]],ts[order[3*len(order)//4]]
        slope = (math.sin(t1-t0))/(g1-g0+1e-9)  # rough
        # straightness: max deviation of theta from linear fit
        print(f"  {name:22s}: {len(pts)} pts; gamma[{round(min(gm))},{round(max(gm))}] "
              f"theta[{round(min(tm))},{round(max(tm))}] deg")


def rotate(words):
    """90-degree tangle rotation (the 'turn' of rational-tangle calculus, s -> -1/s):
    cyclic relabel of the 4 boundary punctures (X1,X2,X3,X4) -> (X2,X3,X4,X1).
    NOTE (2026-07-14): this produces a valid straight-line curve for each rotated
    integer tangle, but the slope reads steep (the rotation swaps the gamma/theta
    roles in pillowcase_coords), so the 1/q vs q labeling still needs pinning
    against the assembled knot (does P(-2,3,5) recover T(3,5)'s 9 earring reps?).
    Convention TBD -- see RESEARCH_LOG sec 21."""
    return [words[1], words[2], words[3], words[0]]
