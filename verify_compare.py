#!/usr/bin/env python3
"""Same-model comparison, not an optimized factory or hardware simulation.
Requires Python >= 3.10 and mpmath.
Exact: syndrome propagation, endpoint comparison, one-fault enumeration.
Numerical: end-to-end expected pooled/retry copy costs at fixed accuracy.
"""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import json
import mpmath as mp

E={1:1,2:2,3:1,4:4,5:1,6:2,7:1}
R3=tuple(("S",i) for i in (1,2,4))
R4=tuple(("S",i) for i in (1,2,4,7))
H12=tuple(zip(("S","S","S","R","S","S","S","S","S","S","S","S"),
              (1,2,4,7,1,3,7,6,3,4,1,6)))

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def output(r,s,ell):
    a=(r & ell).bit_count()%2
    b=(s & ell).bit_count()%2
    return r^s^(E[ell] if b else 0) if a==b else None

def combine(P,Qdist,ell):
    out=[0]*8
    for r in range(8):
        for s in range(8):
            u=output(r,s,ell)
            if u is not None:
                out[u]+=P[r]*Qdist[s]
    return out

def run(F,schedule,lam=None):
    raw=[F]+[(1-F)/7]*7
    P=list(raw); cost=1; Z=1; leaves=1
    for mode,ell in schedule:
        w=combine(P,P if mode=="S" else raw,ell)
        A=sum(w)
        require(A>0,"zero acceptance")
        P=[a/A for a in w]
        cost=(2*cost if mode=="S" else cost+1)/A
        Z=A*(Z**2 if mode=="S" else Z)
        leaves=2*leaves if mode=="S" else leaves+1
        if lam is not None:
            P=[(1-lam)*p+lam/7*(1-p) for p in P]
    return P[0],cost,Z,leaves

def text(x,digits=35):
    if isinstance(x,Q):
        x=mp.mpf(x.numerator)/x.denominator
    return mp.nstr(x,digits)

def factory(Fs,delta_s,round_schedule,prefix=None):
    F=mp.mpf(Fs); cost=mp.mpf(1); rounds=0
    if prefix is not None:
        F,cost,_,_=run(F,prefix)
    while 1-F>mp.mpf(delta_s):
        require(rounds<50,"iteration cap exceeded")
        Fnew,multiplier,_,_=run(F,round_schedule)
        require(Fnew>F,"does not improve")
        F=Fnew; cost*=multiplier; rounds+=1
    return {"rounds_after_prefix":rounds,
            "final_infidelity":text(1-F),
            "expected_raw_copies":text(cost)}

def fault_rows(schedule):
    word=[ell for _,ell in schedule]; rows=[]
    for level in range(len(word)):
        locs=2**(len(word)-level-1)
        wrong=rejected=correct=0
        for loc in range(locs):
            for d in range(1,8):
                P=[0]*(2**len(word))
                for k,ell in enumerate(word):
                    nxt=[]
                    for j in range(0,len(P),2):
                        u=output(P[j],P[j+1],ell)
                        if u is None:
                            break
                        if k==level and j//2==loc:
                            u ^= d
                        nxt.append(u)
                    if len(nxt)!=len(P)//2:
                        P=None
                        break
                    P=nxt
                if P is None: rejected+=1
                elif P[0]: wrong+=1
                else: correct+=1
        rows.append([locs,wrong,rejected,correct])
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=Path("comparison_results.json"),
                    help="destination for the machine-readable results")
    args=ap.parse_args()
    mp.mp.dps=80
    results={"assumptions":{
        "state":"isotropic even-parity independent raw inputs",
        "cost":"expected raw copies with independent retry/pooling at each comparison",
        "not_counted":"storage noise, elapsed time, spatial routing, microscopic gate costs",
        "stopping":"whole prescribed rounds; first output meeting target; not optimized",
        "original_paper_cost_caveat":
        "Bravyi Eq. (33) instead counts acceptance of a whole preassigned round. "
        "We grant elementary-comparison pooling to BOTH recurrences here."}}
    results["low_fidelity"]={}
    for Fs in (".5995",".5997"):
        entry={}
        for name,schedule in (("R3",R3),("H12",H12)):
            f,c,z,n=run(Q(Fs),schedule)
            entry[name]={"Fout":text(f),"cost":text(c),
                         "accepted_tree_leaves":n,"accepted_tree_nodes":n-1}
            if name=="H12":
                require(f>Q(664887507913,1099511627776),"endpoint basin entry failed")
        results["low_fidelity"][Fs]=entry
    results["finite_accuracy_numerical"]={}
    for Fs in (".65",".9"):
        row={}
        for name,sched,pre in (("R3",R3,None),("R4",R4,None),("H12_then_R4",R4,H12)):
            row[name]=factory(Fs,"1e-6",sched,pre)
            with mp.workdps(110):
                check=factory(Fs,"1e-6",sched,pre)
                require(check["rounds_after_prefix"]==row[name]["rounds_after_prefix"],
                        "precision changed stopping rule")
                require(abs(mp.mpf(check["expected_raw_copies"])/
                            mp.mpf(row[name]["expected_raw_copies"])-1)<mp.mpf("1e-30"),
                        "cost unstable across precisions")
        results["finite_accuracy_numerical"][Fs]=row
    # Low-fidelity endpoints of the complete prefix+R4 factory in
    # Eq. (endpointfactorycost). These costs are numerical, not certified
    # optimal factory costs or exact rational enclosures.
    results["endpoint_factories_numerical"] = {}
    for Fs, expected_rounds in ((".5995", 8), (".5997", 5)):
        row = factory(Fs, "1e-6", R4, H12)
        require(row["rounds_after_prefix"] == expected_rounds,
                f"unexpected endpoint factory stopping round for {Fs}")
        with mp.workdps(110):
            again = factory(Fs, "1e-6", R4, H12)
            require(again["rounds_after_prefix"] == expected_rounds,
                    "endpoint factory stopping changed with precision")
            require(abs(mp.mpf(again["expected_raw_copies"])/
                        mp.mpf(row["expected_raw_copies"])-1) < mp.mpf("1e-30"),
                    "endpoint factory cost unstable")
        results["endpoint_factories_numerical"][Fs] = row
    results["single_postcheck_faults"]={}

    for name,sched,expected in (
        ("R3",R3,[[4,4,24,0],[2,6,8,0],[1,7,0,0]]),
        ("R4",R4,[[8,0,56,0],[4,4,24,0],[2,6,8,0],[1,7,0,0]])):
        rows=fault_rows(sched)
        require(rows==expected,"fault counts differ")
        a=sum(r[1] for r in rows)
        deriv=mp.diff(lambda lam:1-run(mp.mpf(1),sched,lam)[0],mp.mpf(0))
        require(a==17 and abs(deriv-mp.mpf(17)/7)<mp.mpf("1e-60"),
                "one-fault coefficient differs")
        results["single_postcheck_faults"][name]={
            "columns":["locations","accepted_wrong","rejected","accepted_correct"],
            "rows":rows,"first_order_infidelity_coefficient":str(Q(a,7))}
    dest=args.output
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(results,indent=2),encoding="utf-8")
    print(json.dumps(results,indent=2))
    print("ALL COMPARISON CHECKS PASSED")

if __name__=="__main__":
    main()
