"""Stage K: complete tie/loss classification of every single mutant of the champion.

Stage G screened all 800x3 mutants forward (champion in slot a) and none beat it.
This adds the reverse query so each mutant gets a three-way verdict, which is what
we need before shipping any of them in a submission.
"""
import json, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge
panel=dict(load_panel())
rank=json.load(open("checkpoints/stage_b_ranking.json"))
CH=rank[0]["name"]; champ=panel[CH]
muts={}
for i in range(800):
    for b in "ACGT":
        if b==champ[i]: continue
        muts[f"M0_{i:03d}{b}"]=champ[:i]+b+champ[i+1:]
print(f"{len(muts)} single mutants", flush=True)
J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(it):
    n,s=it; key=f"{n}|{CH[:36]}"
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(s,champ,n[:36],CH[:36])
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%800==0: print(f"  {done[0]}/{len(muts)} paid={J.paid}",flush=True); J.save()
    return n,v
res=list(ThreadPoolExecutor(max_workers=12).map(run, muts.items()))
J.save()
c=J.cache; from collections import Counter
verd={}
for n,rev in res:
    fwd=c.get(f"{CH[:36]}|{n}")
    if fwd is None: verd[n]="no_forward"; continue
    verd[n]= "BEATS" if (not fwd and rev) else "loses" if (fwd and not rev) else "tie"
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print("verdicts:", dict(Counter(verd.values())))
losers=[n for n,v in verd.items() if v=="loses"]
if losers:
    print(f"\n{len(losers)} mutants STRICTLY WORSE than the champion -- positions:")
    pos=sorted({int(n[3:6]) for n in losers})
    print(" ", pos[:60], "..." if len(pos)>60 else "")
json.dump(verd, open("checkpoints/stage_k_mutant_verdicts.json","w"), indent=1)
