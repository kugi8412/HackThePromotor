"""Stage N: can any champion-tying variant BEAT a plateau member?

The champion only ties with the other 35 plateau sequences. A variant that beats one
from the disadvantaged slot is therefore strictly stronger than the champion on that
comparison -- the only route left to breaking the ceiling. Anchors sit in slot `a`.
"""
import json, threading, time, collections
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge, tata, gc

panel=dict(load_panel()); rank=json.load(open("checkpoints/stage_b_ranking.json"))
yard=json.load(open("checkpoints/stage_i.json"))
tierA=[r["name"] for r in rank if r["wins"]>0]
tierB=sorted([r["name"] for r in rank if r["wins"]==0 and r["losses"]==0], key=lambda n:-yard.get(n,0))
anchors=[tierA[1], tierA[4], tierB[0]]          # two tier-A, one tier-B, all on the ceiling
print("anchors:", anchors)

m=json.load(open("checkpoints/stage_m.json")); V=m["verdicts"]; par=m["parent"]
p=r"C:\Users\leszk\Downloads\output_5changes.fasta"
recs=[]; n=None; buf=[]
for line in open(p):
    line=line.strip()
    if line.startswith(">"):
        if n: recs.append((n,"".join(buf)))
        n,buf=line[1:],[]
    else: buf.append(line)
recs.append((n,"".join(buf)))
orig={panel[x] for x in (tierA+tierB)[:20]}
seen={}
for h,s in recs:
    if s in orig: continue
    seen.setdefault(s,h)
seqs=[(f"V{i:05d}",s) for i,s in enumerate(seen)]
ties=[(nm,s) for nm,s in seqs if V.get(nm)=="tie"]
print(f"champion-tying variants: {len(ties)}   queries: {len(ties)*len(anchors)}", flush=True)

J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(it):
    nm,s=it; wins=[]
    for a in anchors:
        k=f"{a}|{nm}"
        with lock: hit=J.cache.get(k)
        if hit is None:
            v=J.h.sedzia(panel[a],s,a[:36],nm)
            with lock: J.cache[k]=v; J.paid+=1
        else:
            v=hit
            with lock: J.hits+=1
        if not v: wins.append(a)
    with lock:
        done[0]+=1
        if done[0]%1500==0: print(f"  {done[0]}/{len(ties)} paid={J.paid} {done[0]/(time.time()-t0):.0f}/s",flush=True); J.save()
    return nm,wins
res=list(ThreadPoolExecutor(max_workers=12).map(run, ties))
J.save()
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
dist=collections.Counter(len(w) for _,w in res)
print("anchors beaten per variant:", dict(sorted(dist.items())))
win=[(nm,w) for nm,w in res if w]
print(f"\nvariants beating >=1 plateau anchor: {len(win)}")
smap=dict(ties)
if win:
    win.sort(key=lambda x:-len(x[1]))
    print(f"\n{'variant':10s} {'anchors beaten':>14s}  {'parent':16s} {'TATA':>4s} {'GC':>6s}")
    for nm,w in win[:25]:
        s=smap[nm]
        print(f"{nm:10s} {len(w):14d}  {par.get(nm,'?'):16s} {tata(s):4d} {gc(s):6.3f}")
json.dump({"winners":{nm:smap[nm] for nm,w in win},"anchors_beaten":{nm:w for nm,w in win}},
          open("checkpoints/stage_n.json","w"))
