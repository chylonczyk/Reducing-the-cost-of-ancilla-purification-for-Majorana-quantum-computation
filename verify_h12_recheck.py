#!/usr/bin/env python3
"""Independent recheck: unnormalized integer weights, exact state/cost identities,
and a separate probability-space operation-noise calculation.
No assert statements. Numerical roots are not interval certificates.
Run: python verify_h12_recheck.py --output audit_output [--reference PATH]
Dependencies: Python >=3.10, SymPy, mpmath.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, sys
from fractions import Fraction as Q
from dataclasses import dataclass
from pathlib import Path
import mpmath as mp
import sympy as sp

H12 = (("S",1),("S",2),("S",4),("R",7),("S",1),("S",3),
       ("S",7),("S",6),("S",3),("S",4),("S",1),("S",6))
W25 = tuple(("S",i) for i in [6,1,5,1,2,7,4,5,2]+[1,6,5,2]*4)
W14 = tuple(("S",i) for i in [6,1,5,1,2,7,4,5,2,1,6,5,2,1])
R4 = tuple(("S",i) for i in [1,2,4,7])
E = (0,1,2,1,4,1,2,1)
LABELS = tuple((a,b,c) for a in (0,1) for b in (0,1) for c in (0,1))
INDEX = {s:i for i,s in enumerate(LABELS)}
LO = Q(1329775015825,2199023255552)
HI = Q(664887507913,1099511627776)

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def parity(a,b):
    return sum(x*y for x,y in zip(a,b)) % 2

def transitions(ell, e=None):
    require(type(ell) is int and 1 <= ell <= 7, "invalid check")
    e = E[ell] if e is None else e
    require(type(e) is int and 0 <= e < 8 and parity(LABELS[ell],LABELS[e])==1,
            "invalid section")
    out=[]
    for r,rv in enumerate(LABELS):
        for s,sv in enumerate(LABELS):
            a,b=parity(LABELS[ell],rv),parity(LABELS[ell],sv)
            if a==b:
                t=tuple((rv[j]+sv[j]+b*LABELS[e][j])%2 for j in range(3))
                out.append((r,s,INDEX[t]))
    return tuple(out)
TRANS={i:transitions(i) for i in range(1,8)}

def combine(left,right,ell):
    w=[0]*8
    for r,s,t in TRANS[ell]:
        w[t] += left[r]*right[s]
    return w

def fractions_fourier(P):
    return [sum(((-1)**parity(v,s))*p for s,p in zip(LABELS,P))
            for v in LABELS]

def mpq(q):
    return mp.mpf(q.numerator)/q.denominator

def enclose(q,digits):
    d=10**digits
    a=q.numerator*d//q.denominator
    b=-((-q.numerator*d)//q.denominator)
    def render(n):
        sign="-" if n<0 else ""
        z=str(abs(n)).zfill(digits+1)
        return sign+(z[:-digits]+"."+z[-digits:] if digits else z)
    require(Q(a,d)<=q<=Q(b,d),"invalid decimal enclosure")
    return render(a),render(b)

def exact_tree(F,schedule=H12,reference=None):
    require(isinstance(F,Q) and 0<=F<=1,"invalid raw fidelity")
    raw=[7*F.numerator]+[F.denominator-F.numerator]*7
    weights=list(raw); leaves=1; cost=Q(1); stages=[]; accepted=[]
    for j,(mode,ell) in enumerate(schedule,1):
        require(mode in ("S","R"),"invalid mode")
        right=weights if mode=="S" else raw
        before=sum(weights)*sum(right)
        new=combine(weights,right,ell)
        total=sum(new)
        require(total>0,"zero acceptance")
        A=Q(total,before)
        if mode=="S": cost=2*cost/A; leaves*=2
        else: cost=(cost+1)/A; leaves+=1
        weights=new
        P=[Q(v,total) for v in weights]
        require(sum(P)==1 and min(P)>=0,"invalid state")
        accepted.append(A)
        if reference is not None:
            f,c,n,aa,x=reference.run_fourier(F,schedule[:j])
            require(f==P[0] and c==cost and n==leaves and aa==accepted,
                    f"reference mismatch at stage {j}")
            require(list(x)==fractions_fourier(P),f"state mismatch at stage {j}")
        stages.append({"stage":j,"mode":mode,"ell":ell,"F":P[0],"A":A,
                       "C":cost,"N":leaves,"P":P})
    # Independent reverse-demand accounting (number of requested input streams).
    demand=Q(1); direct_raw=Q(0); multiplicities=[]; multiplicity=1
    for j in range(len(schedule)-1,-1,-1):
        mode,ell=schedule[j]
        multiplicities.append(multiplicity)
        attempts=demand/accepted[j]
        if mode=="S": demand=2*attempts; multiplicity*=2
        else: demand=attempts; direct_raw+=attempts
    require(demand+direct_raw==cost,"reverse-demand cost differs")
    require(sum(multiplicities)==leaves-1,"full binary tree node count differs")
    # The full-tree acceptance weight has denominator d^N, not the factory cost.
    Z=Q(sum(weights),sum(raw)**leaves)
    z=Q(1)
    for (mode,_),A in zip(schedule,accepted):
        z=A*(z*z if mode=="S" else z)
    require(z==Z,"full-tree acceptance recursion mismatch")
    Fout=Q(weights[0],sum(weights))
    rawbytes=(Fout.numerator.to_bytes((Fout.numerator.bit_length()+7)//8,"big")
              +b"/"+Fout.denominator.to_bytes((Fout.denominator.bit_length()+7)//8,"big"))
    return Fout,cost,stages,Z,hashlib.sha256(rawbytes).hexdigest()

def recurrence_check():
    # Reconstruct R4 independently as integer leaf-weight polynomials.
    z,F=sp.symbols("z F")
    raw=[sp.Poly(1,z)]+[sp.Poly(z,z)]*7
    cur=raw
    for _,ell in R4: cur=combine(cur,cur,ell)
    def eval_poly(p):
        return sum(int(p.nth(w))*F**(16-w)*((1-F)/7)**w for w in range(17))
    good=eval_poly(cur[0]); bad=sum(eval_poly(p) for p in cur[1:])
    Z=sp.expand(good+bad)
    D=sp.expand(7**16*Z)
    g=sp.div(sp.Poly(sp.expand(7**16*(good-F*Z)),F),
             sp.Poly(-(F-1)*(8*F-1),F))
    require(g[1].is_zero,"fixed point factorization has remainder")
    G=g[0]
    require(G.degree()==15,"wrong fixed-point degree")
    seq=sp.sturm(G.as_expr(),F)
    def variations(x):
        q=sp.Rational(x.numerator,x.denominator) if isinstance(x,Q) else x
        signs=[sp.sign(v.subs(F,q)) for v in seq]
        signs=[s for s in signs if s!=0]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    require((variations(0),variations(1))==(8,7),"Sturm count differs")
    require(G.eval(sp.Rational(LO.numerator,LO.denominator))<0 and
            G.eval(sp.Rational(HI.numerator,HI.denominator))>0,
            "root enclosure failed")
    require(G.eval(1)>0,"wrong upper-side sign")
    # D>0: nonnegative accepted coefficients, all-zero and all-fault terms >0.
    require(all(p.nth(w)>=0 for p in cur for w in range(17)),"negative enumerator")
    require(cur[0].nth(0)==1 and sum(p.nth(16) for p in cur)>0,
            "structural positivity failed")
    require([sum(p.nth(w) for p in cur[1:]) for w in range(4)]==[0,0,0,384],
            "low-weight error coefficients differ")
    eps=sp.symbols("eps")
    series=sp.series((bad/Z).subs(F,1-eps),eps,0,5).removeO()
    require(sp.expand(series-sp.Rational(384,343)*eps**3
                      -sp.Rational(11520,2401)*eps**4)==0,"Taylor law differs")
    return {"Sturm_variations":[8,7],"root_lower":str(LO),"root_upper":str(HI)}

def basis_check():
    # Exhaustive bilinear numerator test, including unequal deterministic inputs.
    for ell in range(1,8):
        for e in range(8):
            if parity(LABELS[ell],LABELS[e])!=1: continue
            mapping={(r,s):u for r,s,u in transitions(ell,e)}
            for r in range(8):
                xr=[(-1)**parity(v,LABELS[r]) for v in LABELS]
                for s in range(8):
                    xs=[(-1)**parity(v,LABELS[s]) for v in LABELS]
                    A=Q(1+xr[ell]*xs[ell],2)
                    require(A==int((r,s) in mapping),"basis acceptance mismatch")
                    for v in range(8):
                        ve=parity(LABELS[v],LABELS[e])
                        if ve==0:
                            num=Q(xr[v]*xs[v]+xr[v^ell]*xs[v^ell],2)
                        else:
                            num=Q(xr[v]*xs[v^ell]+xr[v^ell]*xs[v],2)
                        exp=((-1)**parity(LABELS[v],LABELS[mapping[(r,s)]])) if A else 0
                        require(num==exp,"basis output mismatch")
    return 28*64

def prob_mp(F,lam,schedule):
    raw=[F]+[(1-F)/7]*7
    P=list(raw)
    for mode,ell in schedule:
        new=combine(P,P if mode=="S" else raw,ell)
        A=mp.fsum(new)
        require(A>0,"numerical zero acceptance")
        P=[(1-lam)*v/A+(lam/7)*(1-v/A) for v in new]
    return P[0]

def bisect(fun,a,b,steps=150):
    a,b=mp.mpf(a),mp.mpf(b)
    fa,fb=fun(a),fun(b)
    require(fa*fb<0,"numerical root not bracketed")
    for _ in range(steps):
        c=(a+b)/2; fc=fun(c)
        if fc==0: return c
        if fa*fc<0: b=c
        else: a=c; fa=fc
    return (a+b)/2

def noise_recheck():
    mp.mp.dps=80
    # Solve simultaneous equations with probability-space updates, then
    # independently bracket the local roots and check derivatives/residuals.
    fstar=bisect(lambda F:prob_mp(F,mp.mpf(0),R4)-F,".60",".61",200)
    fcross=bisect(lambda F:prob_mp(F,mp.mpf(0),H12)-fstar,
                  ".59949",".5995",200)
    results={}
    specs=[("H12_5997",H12,".5997","0.000129318"),
           ("H12_5995",H12,".5995","0.000005655"),
           ("W25_5997",W25,".5997","0.0000001007"),
           ("W14_6000",W14,".6","0.0000251686")]
    for name,sched,Fs,guess in specs:
        F=mp.mpf(Fs)
        def equations(b,lam):
            return (prob_mp(b,lam,R4)-b,prob_mp(F,lam,sched)-b)
        b,lam=mp.findroot(equations,(mp.mpf(".605"),mp.mpf(guess)),tol=mp.mpf("1e-65"))
        require(0<lam<mp.mpf(".001"),"unexpected local crossing")
        require(max(abs(v) for v in equations(b,lam))<mp.mpf("1e-60"),
                "large root residual")
        d=mp.diff(lambda z:prob_mp(z,lam,R4),b)
        require(d>1,"basin boundary not unstable")
        upper=mp.findroot(lambda z:prob_mp(z,lam,R4)-z,(mp.mpf(".99"),mp.mpf(1)))
        require(abs(mp.diff(lambda z:prob_mp(z,lam,R4),upper))<1,
                "upper fixed point not stable")
        # Check smaller and larger lambda sides at 90 and 110 percent of crossing.
        margins=[]
        for factor in (mp.mpf(".9"),mp.mpf("1.1")):
            ll=factor*lam
            lower=bisect(lambda z:prob_mp(z,ll,R4)-z,".60",".62",120)
            margins.append(prob_mp(F,ll,sched)-lower)
        require(margins[0]>0 and margins[1]<0,"crossing direction differs")
        # Repeat joint solve with 110 digits.
        with mp.workdps(110):
            b2,l2=mp.findroot(equations,(b,lam),tol=mp.mpf("1e-90"))
            require(abs(l2-lam)<mp.mpf("1e-60"),"precision instability")
        results[name]={"Fraw":Fs,"lambda":mp.nstr(lam,55),"lower_fixed_point":mp.nstr(b,45),
                       "R4_slope":mp.nstr(d,30),"signs_at_90_110_percent":["+","-"]}
        print(name,"local crossing =",mp.nstr(lam,30),flush=True)
    results["H12_noiseless_crossing"]=mp.nstr(fcross,55)
    results["lambda_ratio_same_F"]=mp.nstr(mp.mpf(results["H12_5997"]["lambda"])/
                                        mp.mpf(results["W25_5997"]["lambda"]),40)
    # Separate structural one-final-channel bound for 0<=lambda<=7/8.
    for lam in (Q(0),Q(1,10000),Q(1,10),Q(7,8)):
        for f in (Q(0),Q(1,2),Q(1)):
            require((1-Q(8,7)*lam)*f+lam/7 <= 1-lam,"last-channel bound fails")
    return results


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        require(self.lo <= self.hi, "reversed interval")

    def __add__(self, other):
        if not isinstance(other, Interval): other=Interval(Q(other),Q(other))
        return Interval(self.lo+other.lo,self.hi+other.hi)

    __radd__=__add__

    def __mul__(self, other):
        if not isinstance(other, Interval): other=Interval(Q(other),Q(other))
        ends=(self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi)
        return Interval(min(ends),max(ends))

    __rmul__=__mul__

    def __truediv__(self, other):
        if not isinstance(other, Interval): other=Interval(Q(other),Q(other))
        require(other.lo>0,"interval denominator not strictly positive")
        return self*Interval(1/other.hi,1/other.lo)

    def rounded(self, digits=80):
        d=10**digits
        a=self.lo.numerator*d//self.lo.denominator
        b=-((-self.hi.numerator*d)//self.hi.denominator)
        return Interval(Q(a,d),Q(b,d))

def interval_baseline():
    F=Q("0.5997")
    raw=[Interval(F,F)]+[Interval((1-F)/7,(1-F)/7)]*7
    P=raw; cost=Interval(Q(1),Q(1))
    for _,ell in W25:
        new=combine(P,P,ell)
        A=sum(new)
        require(A.lo>0,"W25 interval acceptance failed")
        P=[(v/A).rounded() for v in new]
        cost=(cost*2/A).rounded()
    y=Interval(Q(1),Q(1))/cost
    require(P[0].lo > HI,"independent W25 interval misses basin")
    require(Q("0.60550725738718132769784")<P[0].lo and
            P[0].hi<Q("0.60550725738718132769785"),
            "independent W25 fidelity enclosure differs")
    require(Q("0.00000000000005684665")<y.lo and
            y.hi<Q("0.00000000000005684666"),
            "independent W25 yield enclosure differs")
    return P[0],cost,y

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="audit_output")
    ap.add_argument("--reference",type=Path)
    ap.add_argument("--skip-noise",action="store_true")
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    ref=None
    if args.reference:
        spec=importlib.util.spec_from_file_location("reference_h12",args.reference)
        ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    mp.mp.dps=80
    results={"basis_pairs_verified":basis_check(),"R4":recurrence_check(),"exact":{}}
    for Fs in ("0.5995","0.5997"):
        F=Q(Fs)
        f,c,rows,Z,digest=exact_tree(F,H12,ref)
        require(f>HI,"H12 does not enter basin")
        if Fs=="0.5995":
            require((f.numerator.bit_length(),f.denominator.bit_length())==(25520,25521),
                    "fraction lengths differ")
            require(digest=="2ccfebe48dd4c020d42e3991fb9e354c9e9a1e692f17b0bda28aed044ca43a51",
                    "regression fingerprint differs")
            require(16*(1-F)/(15-8*F)==Q(1602,2551),"depolarization conversion differs")
        result={"Fout":enclose(f,30),"cost":enclose(c,9),"yield":enclose(1/c,25),
                "leaves":rows[-1]["N"],"nodes":rows[-1]["N"]-1,
                "margin_above_Fstar_upper":enclose(f-HI,25),"sha256":digest,
                "Fout_numerator_hex":format(f.numerator,"x"),
                "Fout_denominator_hex":format(f.denominator,"x")}
        results["exact"][Fs]=result
        with (out/("stages_"+Fs.replace(".","")+".csv")).open("w",newline="") as fp:
            writer=csv.writer(fp)
            writer.writerow(["stage","mode","ell","F_low","F_high","A_low","A_high",
                             "C_low","C_high","raw_leaves"])
            for r in rows:
                writer.writerow([r["stage"],r["mode"],r["ell"],*enclose(r["F"],18),
                                 *enclose(r["A"],18),*enclose(r["C"],9),r["N"]])
        print("H12",Fs,":",result["Fout"],"cost",result["cost"],flush=True)
        # Explicit nonrecursive cost formula for this one-injection schedule.
        AA=[r["A"] for r in rows]
        closed=(Q(2048)/math.prod(AA[:3])+256)/math.prod(AA[3:])
        require(closed==c,"closed-form H12 cost differs")
    # Independently re-enclose the W25 benchmark in exact probability-space
    # interval arithmetic; no stored decimal result is used to propagate it.
    fw,cw,yw=interval_baseline()
    results["W25_independent_interval"]={
        "Fout":[str(fw.lo),str(fw.hi)],"cost":[str(cw.lo),str(cw.hi)],
        "yield":[str(yw.lo),str(yw.hi)]}
    h_yield=1/exact_tree(Q("0.5997"))[1]
    results["yield_ratio_from_recomputed_W25_interval"]={
        "lower":enclose(h_yield/yw.hi,15)[0],"upper":enclose(h_yield/yw.lo,15)[1]}
    print("Same-F yield ratio:",results["yield_ratio_from_recomputed_W25_interval"],
          flush=True)
    if not args.skip_noise: results["noise_numerical"]=noise_recheck()
    results["environment"]={"python":sys.version,"sympy":sp.__version__,"mpmath":mp.__version__}
    (out/"recheck_results.json").write_text(json.dumps(results,indent=2))
    print("ALL INDEPENDENT RECHECKS PASSED",flush=True)

if __name__=="__main__":
    main()
