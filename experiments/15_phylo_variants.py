"""Stage M: screen the supplied 5-point-mutation variants against the champion.

Both orders, because the question is "ties OR wins" and a single query from the
champion's advantaged slot cannot separate a tie from a loss. Parent assignment is
by minimum Hamming distance to the 20 ranked naturals, since the file's headers
don't record it.
"""
import json, threading, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from judge_scan import load_panel, CachedJudge, tata, gc

panel=dict(load_panel()); rank=json.load(open("checkpoints/stage_b_ranking.json"))
CH=rank[0]["name"]; champ=panel[CH]
tierA=[r["name"] for r in rank if r["wins"]>0]
yard=json.load(open("checkpoints/stage_i.json"))
tierB=sorted([r["name"] for r in rank if r["wins"]==0 and r["losses"]==0], key=lambda n:-yard.get(n,0))
parents=(tierA+tierB)[:20]

p=r"C:\Users\leszk\Downloads\output_5changes.fasta"
recs=[]; n=None; buf=[]
for line in open(p):
    line=line.strip()
    if line.startswith(">"):
        if n: recs.append((n,"".join(buf)))
        n,buf=line[1:],[]
    else: buf.append(line)
recs.append((n,"".join(buf)))
orig={panel[x] for x in parents}
uniq={}
for h,s in recs:
    if s in orig: continue                       # skip the unmodified parents
    uniq.setdefault(s, h)
seqs=[(f"V{i:05d}", s) for i,(s,h) in enumerate(uniq.items())]
print(f"unique variants to screen: {len(seqs)}", flush=True)

A=np.frombuffer("".join(s for _,s in seqs).encode(), dtype=np.uint8).reshape(len(seqs),800)
P=np.frombuffer("".join(panel[x] for x in parents).encode(), dtype=np.uint8).reshape(len(parents),800)
dist=np.stack([(A!=P[k]).sum(1) for k in range(len(parents))],1)
pidx=dist.argmin(1); pdist=dist.min(1)
print(f"distance to assigned parent: min {pdist.min()} median {int(np.median(pdist))} max {pdist.max()}")
import collections
print("variants per parent rank:", dict(sorted(collections.Counter(pidx.tolist()).items())))

J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def ask(a,sa,b,sb):
    k=f"{a}|{b}"
    with lock: hit=J.cache.get(k)
    if hit is not None:
        with lock: J.hits+=1
        return hit
    v=J.h.sedzia(sa,sb,a[:36],b[:36])
    with lock: J.cache[k]=v; J.paid+=1
    return v
def run(it):
    i,(nm,s)=it
    fwd=ask(CH[:36],champ,nm,s)
    rev=ask(nm,s,CH[:36],champ)
    with lock:
        done[0]+=1
        if done[0]%2000==0: print(f"  {done[0]}/{len(seqs)} paid={J.paid} {done[0]/(time.time()-t0):.0f}/s",flush=True); J.save()
    return nm,("BEATS" if (not fwd and rev) else "loses" if (fwd and not rev) else "tie")
res=list(ThreadPoolExecutor(max_workers=12).map(run, enumerate(seqs)))
J.save()
V=dict(res)
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print("verdict vs champion:", dict(collections.Counter(V.values())))
by=collections.defaultdict(collections.Counter)
for i,(nm,s) in enumerate(seqs): by[parents[pidx[i]]][V[nm]]+=1
print(f"\n{'parent':16s} {'rank':>4s} {'n':>5s} {'BEATS':>6s} {'tie':>6s} {'loses':>6s}")
for k,pn in enumerate(parents):
    c=by[pn]; tot=sum(c.values())
    if tot: print(f"{pn:16s} {k+1:4d} {tot:5d} {c['BEATS']:6d} {c['tie']:6d} {c['loses']:6d}")
win=[nm for nm,v in V.items() if v=="BEATS"]
print(f"\nWINNERS: {len(win)}")
seqmap=dict(seqs)
json.dump({"verdicts":V,"winners":{w:seqmap[w] for w in win},
           "parent":{seqs[i][0]:parents[pidx[i]] for i in range(len(seqs))}},
          open("checkpoints/stage_m.json","w"))
