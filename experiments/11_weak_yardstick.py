"""Stage I: rank the 36-member plateau against a WEAKER yardstick.

Within-plateau comparisons are all ties, so head-to-head carries no information.
Instead score each plateau member by how many stage-A rejects it can beat FROM THE
DISADVANTAGED SLOT -- reject in slot `a`, plateau member in slot `b`. Putting the
strong sequence in slot `a` would saturate; this way the count spreads.
"""
import json, random, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge

panel=dict(load_panel())
rank=json.load(open("checkpoints/stage_b_ranking.json"))
plateau=[r["name"] for r in rank if r["losses"]==0]
tierA={r["name"] for r in rank if r["wins"]>0}
survivors=set(json.load(open("checkpoints/stage_a_winners.json")))
rejects=[n for n in panel if n not in survivors]
rng=random.Random(7); yard=rng.sample(rejects,100)
jobs=[(p,y) for p in plateau for y in yard]
print(f"{len(plateau)} plateau x {len(yard)} yardstick = {len(jobs)} calls", flush=True)

J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(j):
    p,y=j; key=f"{y}|{p}"                    # yardstick in slot a, plateau in slot b
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(panel[y],panel[p],y,p)
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%1000==0: print(f"  {done[0]}/{len(jobs)} paid={J.paid}",flush=True); J.save()
    return p,(not v)                          # True = plateau member won from slot b
res=list(ThreadPoolExecutor(max_workers=12).map(run,jobs))
J.save()
from collections import defaultdict
sc=defaultdict(int)
for p,w in res: sc[p]+=w
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s\n")
old={r["name"]:r["wins"] for r in rank}
print(f"{'sequence':16s} {'tier':4s} {'old W/115':>9s} {'new W/100':>9s}")
for p in sorted(plateau,key=lambda x:-sc[x]):
    print(f"{p:16s} {'A' if p in tierA else 'B':4s} {old[p]:9d} {sc[p]:9d}")
vals=[sc[p] for p in plateau]
print(f"\nspread: min {min(vals)}  max {max(vals)}  distinct values {len(set(vals))}/36")
b=[sc[p] for p in plateau if p not in tierA]
print(f"tier B spread: min {min(b)} max {max(b)} distinct {len(set(b))}/{len(b)}")
json.dump({p:sc[p] for p in plateau}, open("checkpoints/stage_i.json","w"), indent=1)
