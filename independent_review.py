#!/usr/bin/env python3
"""Independent review of the Majorana-ancilla purification calculation.

Basis: Fojcik's supplied research notes (July 2026) and verification scripts.
This file does NOT call their Fourier-update implementation.
Core checks use exact integer/rational probability-space recursions.

Run:
    python independent_review.py
    python -O independent_review.py
    python independent_review.py --cross-check verify_modelB_gate_noise_reviewed.py

The full 16-Majorana braid network is NOT simulated here: the published
Bravyi syndrome instrument is the explicit physical input.
The state-dependent fault counterexample tests a missing assumption in the
original worst-case statement; it does not invalidate the ideal threshold.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from collections import Counter
from fractions import Fraction as Q
from itertools import product, permutations
from pathlib import Path
import platform
import sys
import sympy as sp

E = {1: 1, 2: 2, 3: 1, 4: 4, 5: 1, 6: 2, 7: 1}
LABELS = (5, 7, 3, 1, 4, 6, 2, 0)
R3 = (1, 2, 4)
R4 = (1, 2, 4, 7)
W25 = (6, 1, 5, 1, 2, 7, 4, 5, 2) + (1, 6, 5, 2) * 4
G7 = (1458176, -884736, 548864, -368640, 138496, -67808, 5440, -6249)
G15 = (8505108987904, -16854525411328, 28467580108800,
       -27066279919616, 18851022503936, -10509663338496,
       4884383924224, -2048766509056, 716242124800, -238242988032,
       58392420352, -20373018624, 4691993600, -2122800736,
       263381416, -149949201)
FSTAR_UPPER = Q(664887507913, 1099511627776)
Matrix = tuple[tuple[int, ...], ...]
Monomial = tuple[int, int]  # real sign, support bitmask
Interval = tuple[Q, Q]
I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)

def dot(a: int, b: int) -> int:
    return (a & b).bit_count() % 2

def pair(r: int, s: int, ell: int, e: int | None = None) -> int | None:
    if e is None:
        e = E[ell]
    if dot(ell, r) != dot(ell, s):
        return None
    return r ^ s ^ (e if dot(ell, s) else 0)

def matvec(m: Matrix, x: int) -> int:
    v = ((x >> 2) & 1, (x >> 1) & 1, x & 1)
    out = [sum(m[i][j] * v[j] for j in range(3)) % 2 for i in range(3)]
    return 4*out[0] + 2*out[1] + out[2]

def matmul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) % 2
                       for j in range(3)) for i in range(3))

def transpose(a: Matrix) -> Matrix:
    return tuple(tuple(a[j][i] for j in range(3)) for i in range(3))

def all_gl() -> list[Matrix]:
    mats = []
    for bits in product((0, 1), repeat=9):
        m = tuple(tuple(bits[3*i:3*i+3]) for i in range(3))
        if len({matvec(m,x) for x in range(8)}) == 8:
            mats.append(m)
    require(len(mats) == 168, "wrong GL(3,2) order")
    return mats

def cmul(a: Monomial, b: Monomial) -> Monomial:
    sa, ma = a
    sb, mb = b
    inv = sum((mb & ((1 << i)-1)).bit_count()
              for i in range(8) if ma & (1 << i))
    return sa * sb * (-1)**inv, ma ^ mb

def cmon(indices: tuple[int, ...], sign: int = 1) -> Monomial:
    out = (sign, 0)
    for p in indices:
        out = cmul(out, (1, 1 << (p-1)))
    return out

def permute_mon(a: Monomial, perm: list[int]) -> Monomial:
    out = (a[0], 0)
    for i in range(8):
        if a[1] & (1 << i):
            out = cmul(out, (1, 1 << perm[i]))
    return out

def even_equivalent(a: Monomial, b: Monomial) -> bool:
    return a == b or a == cmul(b, (1, 255))

def verify_operator_frames() -> dict[str, object]:
    sgens = [cmon((1,2,5,6),-1), cmon((2,3,6,7),-1), cmon((1,2,3,4),-1)]
    derived = []
    for p in range(8):
        bits = [int(cmul((1,1<<p), s) != cmul(s,(1,1<<p))) for s in sgens]
        derived.append(4*bits[0]+2*bits[1]+bits[2])
    require(tuple(derived) == LABELS, "operator label derivation mismatch")
    stabilizers = []
    for v in range(8):
        out = (1,0)
        for i, s in enumerate(sgens):
            if (v >> (2-i)) & 1:
                out = cmul(out,s)
        stabilizers.append(out)
    mats = all_gl()
    offsets = Counter()
    inverses = {m: next(n for n in mats if matmul(m,n) == I3) for m in mats}
    for m in mats:
        perm = [LABELS.index(matvec(m,x)) for x in LABELS]
        inv = sum(perm[i] > perm[j] for i in range(8) for j in range(i+1,8))
        require(inv % 2 == 0, "a GL label permutation is odd")
        mit = transpose(inverses[m])
        signs = []
        for v in range(8):
            lhs = permute_mon(stabilizers[v],perm)
            rhs = stabilizers[matvec(mit,v)]
            if even_equivalent(lhs,rhs):
                signs.append(0)
            else:
                require(even_equivalent(lhs,(-rhs[0],rhs[1])),
                        "stabilizer support/covariance mismatch")
                signs.append(1)
        b = 4*signs[4]+2*signs[2]+signs[1]
        require(all(signs[v] == dot(b,v) for v in range(8)),
                "operator sign offset is not a character")
        offsets[b] += 1
    require(all(offsets[b] == 21 for b in range(8)), "unexpected offset counts")
    # Bravyi's three syndrome shears generate the entire linear group.
    generators = (((1,1,0),(0,1,0),(0,0,1)),
                  ((1,0,0),(0,1,1),(0,0,1)),
                  ((1,0,0),(0,1,0),(1,0,1)))
    seen = {I3}
    todo = [I3]
    while todo:
        m = todo.pop()
        for g in generators:
            n = matmul(g,m)
            if n not in seen:
                seen.add(n)
                todo.append(n)
    require(seen == set(mats), "Bravyi shears do not generate GL")
    planes = {frozenset(x for x in range(8) if dot(v,x) == b)
              for v in range(1,8) for b in (0,1)}
    affine = {tuple(matvec(m,x)^t for x in range(8)) for m in mats for t in range(8)}
    normalizer = {p for p in permutations(range(8))
                  if all(frozenset(p[x] for x in plane) in planes for plane in planes)}
    require(normalizer == affine and len(affine) == 1344, "normalizer mismatch")
    # Derive feed-forward bit changes, not just their Hamming consequences.
    require(LABELS[1] ^ LABELS[2] == 4, "D23 does not flip s1")
    require(LABELS[0] ^ LABELS[1] == 2, "D12 does not flip s2")
    return {"GL_order":168, "affine_order":1344, "offset_histogram":dict(offsets)}

def verify_fourier_basis() -> None:
    for ell in range(1,8):
        for e in range(1,8):
            if dot(ell,e) != 1:
                continue
            for r in range(8):
                for s in range(8):
                    u = pair(r,s,ell,e)
                    acc = int(u is not None)
                    xr = [(-1)**dot(v,r) for v in range(8)]
                    xs = [(-1)**dot(v,s) for v in range(8)]
                    require(2*acc == 1+xr[ell]*xs[ell], "basis acceptance identity")
                    for v in range(8):
                        direct = 0 if u is None else 2*(-1)**dot(v,u)
                        num = (xr[v]*xs[v]+xr[v^ell]*xs[v^ell]
                               if dot(v,e) == 0 else
                               xr[v]*xs[v^ell]+xr[v^ell]*xs[v])
                        require(direct == num, "unnormalized character identity")

def convolution(a: list[int], b: list[int]) -> list[int]:
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out

def count_polynomials(word: tuple[int,...]) -> list[list[int]]:
    polynomials = [[1,0]]+[[0,1] for _ in range(7)]
    for ell in word:
        out = [[0]*(2*len(polynomials[0])-1) for _ in range(8)]
        for r in range(8):
            for s in range(8):
                u = pair(r,s,ell)
                if u is not None:
                    c = convolution(polynomials[r],polynomials[s])
                    out[u] = [a+b for a,b in zip(out[u],c)]
        polynomials = out
    return polynomials

def variations(sequence: list[sp.Expr], x: sp.Symbol, endpoint: int) -> int:
    signs = [sp.sign(p.subs(x,endpoint)) for p in sequence]
    signs = [s for s in signs if s != 0]
    return sum(a != b for a,b in zip(signs,signs[1:]))

def verify_recurrence(word: tuple[int,...], coeffs: tuple[int,...]):
    counts = count_polynomials(word)
    n = 2**len(word)
    correct = counts[0]
    wrong = [sum(c[w] for c in counts[1:]) for w in range(n+1)]
    f = sp.Symbol("F")
    q = (1-f)/7
    g = sum(c*f**(n-w)*q**w for w,c in enumerate(correct))
    b = sum(c*f**(n-w)*q**w for w,c in enumerate(wrong))
    z = sp.expand(g+b)
    recurrence = sp.cancel(g/z)
    factor = sp.Poly.from_list(coeffs,f).as_expr()
    require(sp.cancel(recurrence-f+(f-1)*(8*f-1)*factor/(7**n*z)) == 0,
            "probability-space recurrence factorization mismatch")
    sturm = sp.sturm(factor,f)
    vv = (variations(sturm,f,0),variations(sturm,f,1))
    require(vv[0]-vv[1] == 1, "wrong exact Sturm root count")
    lower,upper = ((Q("0.61560748080966942"),Q("0.61560748080966943"))
                   if n == 8 else
                   (Q(1329775015825,2199023255552),FSTAR_UPPER))
    require(factor.subs(f,lower)*factor.subs(f,upper) < 0, "root bracket fails")
    ep = sp.Symbol("eps")
    order = 4 if n == 8 else 6
    series = sp.series(1-recurrence.subs(f,1-ep),ep,0,order).removeO()
    expected = (sp.Rational(48,49)*ep**2+sp.Rational(960,343)*ep**3 if n == 8
                else sp.Rational(384,343)*ep**3+sp.Rational(11520,2401)*ep**4
                +sp.Rational(218112,16807)*ep**5)
    require(sp.expand(series-expected) == 0, "local error series fails")
    require(recurrence.subs(f,sp.Rational(3,4)) > sp.Rational(3,4),
            "exact gain sign fails")
    print(f"[ok] R{len(word)} direct enumerators, factorization, Sturm {vv}, Taylor")
    return counts, f, recurrence, z

def add(a: Interval,b: Interval) -> Interval:
    return a[0]+b[0],a[1]+b[1]

def mul(a: Interval,b: Interval) -> Interval:
    vals = [x*y for x in a for y in b]
    return min(vals),max(vals)

def divide(a: Interval,b: Interval) -> Interval:
    require(b[0] > 0, "interval denominator must be strictly positive")
    vals = [x/y for x in a for y in b]
    return min(vals),max(vals)

def rounded(a: Interval,grid: int) -> Interval:
    require(grid > 0, "grid denominator must be positive")
    lo,hi = a[0]*grid,a[1]*grid
    return (Q(lo.numerator//lo.denominator,grid),
            Q(-((-hi.numerator)//hi.denominator),grid))

def decimal_bound(q: Q,digits: int,up: bool=False) -> str:
    """Directed display for arbitrary signs; no floating-point conversion."""
    require(digits >= 0, "negative decimal count")
    scale = 10**digits
    y = q*scale
    n = -((-y.numerator)//y.denominator) if up else y.numerator//y.denominator
    sign = "-" if n < 0 else ""
    n = abs(n)
    if digits == 0:
        return sign+str(n)
    return f"{sign}{n//scale}.{n%scale:0{digits}d}"

def probability_interval_word() -> tuple[Interval,Interval]:
    grid = 10**100
    f = Q(5997,10000)
    probs = [(f,f)]+[((1-f)/7,(1-f)/7) for _ in range(7)]
    acceptance_product = (Q(1),Q(1))
    for ell in W25:
        raw = [(Q(0),Q(0)) for _ in range(8)]
        for r in range(8):
            for s in range(8):
                u = pair(r,s,ell)
                if u is not None:
                    raw[u] = add(raw[u],mul(probs[r],probs[s]))
        acc = (Q(0),Q(0))
        for p in raw:
            acc = add(acc,p)
        require(acc[0] > 0, "zero acceptance interval")
        acceptance_product = rounded(mul(acceptance_product,acc),grid)
        probs = [rounded(divide(p,acc),grid) for p in raw]
    fidelity = probs[0]
    yield_interval = divide(acceptance_product,(Q(2**25),Q(2**25)))
    require(fidelity[0] > FSTAR_UPPER, "direct certificate fails basin entry")
    require(Q("0.60550725738718132769784") < fidelity[0]
            <= fidelity[1] < Q("0.60550725738718132769785"),
            "direct certificate outside compact claimed interval")
    print("[ok] W25 direct probability-space rational interval certificate")
    print("fidelity lower =",decimal_bound(fidelity[0],85))
    print("fidelity upper =",decimal_bound(fidelity[1],85,True))
    print("yield lower =",decimal_bound(yield_interval[0],20))
    print("yield upper =",decimal_bound(yield_interval[1],20,True))
    return fidelity,yield_interval

def check_review_findings(f,recurrence,z) -> dict[str,str]:
    ff = Q(3,4)
    eta = Q(1,10**6)
    zz = Q(z.subs(f,sp.Rational(3,4)))
    rr = Q(recurrence.subs(f,sp.Rational(3,4)))
    require(eta < ff**16, "counterexample fault mass unavailable")
    actual = (zz*rr-eta)/zz
    claimed = (1-eta)*zz*rr/((1-eta)*zz+eta)
    require(actual < claimed, "state-dependent counterexample fails")
    print("[ok] Missing ideal-branch assumption: exact counterexample established")
    print("actual fidelity =",decimal_bound(actual,20))
    print("invalid claimed lower =",decimal_bound(claimed,20))
    # Four raw blocks, two uniform checks: finite batch floors change the yield.
    finite_batch_per_raw = Q(1,2)**2 * Q(1,2) / 4
    asymptotic_pool = Q(1,4)*Q(1,2)**2
    require(finite_batch_per_raw == Q(1,32) < asymptotic_pool == Q(1,16),
            "finite batch example incorrect")
    x = sp.Symbol("x")
    tk = (2*x+12*x*x)/(7*(1+x*x))
    require(tk.subs(x,0) == 0, "missing zero fixed-point check")
    require(set(sp.solve(sp.together(tk-x),x)) == {0,sp.Rational(5,7),1},
            "one-check fixed-point set")
    print("[ok] Finite-batch yield distinction and missing x=0 fixed point")
    return {"actual_counterexample_fidelity":str(actual),
            "invalid_lower_bound":str(claimed),
            "finite_batch_yield":"1/32", "asymptotic_pool_yield":"1/16"}

def verify_target_and_gaussian_closure_example() -> None:
    # Exact Jordan-Wigner matrices for the four-mode target (16-dimensional).
    I=sp.eye(2); X=sp.Matrix([[0,1],[1,0]])
    Y=sp.Matrix([[0,-sp.I],[sp.I,0]]); Z=sp.diag(1,-1)
    cs=[]
    for j in range(4):
        for op in (X,Y):
            cs.append(sp.kronecker_product(*([Z]*j+[op]+[I]*(3-j))))
    def prodmat(ids):
        out=sp.eye(16)
        for i in ids: out=out*cs[i-1]
        return out
    s1=-prodmat((1,2,5,6));s2=-prodmat((2,3,6,7));s3=-prodmat((1,2,3,4))
    parity=prodmat(tuple(range(1,9)))
    p=(sp.eye(16)+parity)*(sp.eye(16)+s1)*(sp.eye(16)+s2)*(sp.eye(16)+s3)/16
    state=sp.zeros(16,1);state[0]=1/sp.sqrt(2);state[15]=1/sp.sqrt(2)
    require(p == state*state.T, "target stabilizer/GHZ convention mismatch")
    # Explicit two-mode Gaussian preparation: K=c1*c3, K^2=-I,
    # exp(pi*K/4)=(I+K)/sqrt(2) maps vacuum to the real BCS Bell state.
    c1=sp.kronecker_product(X,I); c3=sp.kronecker_product(Z,X)
    K=c1*c3
    require(K*K == -sp.eye(4), "quadratic generator square")
    U=(sp.eye(4)+K)/sp.sqrt(2)
    vac=sp.Matrix([1,0,0,0])
    psi=U*vac
    require(psi == sp.Matrix([1,0,0,1])/sp.sqrt(2), "Gaussian BCS preparation")
    rho=psi*psi.conjugate().T
    reduced=sp.Matrix(2,2,lambda i,j: sum(rho[2*i+k,2*j+k] for k in range(2)))
    require(reduced == sp.eye(2)/2 and sp.trace(reduced**2) == sp.Rational(1,2),
            "Gaussian partial-trace example must be mixed")
    print("[ok] Exact target convention; Gaussian partial trace need not stay pure")

def cross_check(module_path: Path, fidelity: Interval, yield_interval: Interval) -> None:
    spec=importlib.util.spec_from_file_location("reviewed_external_verifier",module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load supplied verification module")
    mod=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=mod
    spec.loader.exec_module(mod)
    f,y=mod.verify_w25()
    # Extended/reviewed verifier returns intervals for both quantities.
    require(f.lo <= fidelity[0] <= fidelity[1] <= f.hi,
            "direct probability enclosure is not inside Fourier enclosure")
    require(y.lo <= yield_interval[0] <= yield_interval[1] <= y.hi,
            "direct yield enclosure is not inside Fourier enclosure")
    print("[ok] Independent probability enclosure lies inside Fourier enclosure")

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cross-check",type=Path)
    parser.add_argument("--json",type=Path,help="write exact results as JSON")
    args=parser.parse_args()
    try:
        print(f"[env] Python {platform.python_version()}; SymPy {sp.__version__}")
        frames=verify_operator_frames()
        print("[ok] Operator signs, affine normalizer, and Bravyi-generator closure")
        verify_target_and_gaussian_closure_example()
        verify_fourier_basis()
        print("[ok] 28 checks and 14336 unnormalized character identities")
        c3,f3,r3,z3=verify_recurrence(R3,G7)
        c4,f4,r4,z4=verify_recurrence(R4,G15)
        require([c[2] for c in c3] == [20,16,16,0,16,0,0,0], "R3 weight2 vector")
        require([c[3] for c in c3] == [16,32,32,64,32,64,64,0], "R3 weight3 vector")
        fidelity,yy=probability_interval_word()
        findings=check_review_findings(f4,r4,z4)
        if args.cross_check:
            cross_check(args.cross_check,fidelity,yy)
        if args.json:
            result={"frames":frames,"R3_counts":c3,"R4_counts":c4,
                    "fidelity_interval":[str(x) for x in fidelity],
                    "yield_interval":[str(x) for x in yy],"review_findings":findings}
            args.json.write_text(json.dumps(result,indent=2)+"\n")
        print("ALL INDEPENDENT CHECKS PASSED")
        return 0
    except Exception as exc:
        print(f"INDEPENDENT VERIFICATION FAILED: {exc}",file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
