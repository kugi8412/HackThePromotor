"""Stage F: fine-grained mosaics of the top-10 naturals, 10 bp blocks.

Two arms:
  POS  -- each of the 80 block slots draws independently from the 10 parents, so
          every fragment stays at its native coordinate.
  FREE -- 80 fragments drawn from the pooled 800 (10 parents x 80 slots) and laid
          down in random order, so native position is destroyed.

FREE is the control on POS: if the wagaP positional sensitivity is real, scrambling
coordinates should be markedly worse.
"""
import json, random, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge, tata, gc

panel=dict(load_panel())
rank=json.load(open("checkpoints/stage_b_ranking.json"))
top10=[r["name"] for r in rank[:10]]
CHAMP=top10[0]; champ=panel[CHAMP]
print("parents:", top10, flush=True)
B=10; NB=800//B
frag={p:[panel[p][i*B:(i+1)*B] for i in range(NB)] for p in top10}
pool=[frag[p][i] for p in top10 for i in range(NB)]

rng=random.Random(20260829)
N=3000
seqs={}
for j in range(N):
    combo=[rng.randrange(10) for _ in range(NB)]
    seqs[f"POS_{j:04d}"]="".join(frag[top10[c]][i] for i,c in enumerate(combo))
for j in range(N):
    seqs[f"FREE_{j:04d}"]="".join(rng.sample(pool, NB))
for p in top10[1:]:
    seqs["CTRL_"+p[-6:]]=panel[p]
seqs={k:v for k,v in seqs.items() if len(v)==800}
print(f"{len(seqs)} sequences to screen", flush=True)

J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(it):
    n,s=it; key=f"{CHAMP}|{n}"
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(champ, s, CHAMP[:36], n[:36])
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%1000==0: print(f"  {done[0]}/{len(seqs)} paid={J.paid} {done[0]/(time.time()-t0):.0f}/s",flush=True); J.save()
    return n,(not v)
res=list(ThreadPoolExecutor(max_workers=12).map(run, seqs.items()))
J.save()
from collections import defaultdict
arm=defaultdict(lambda:[0,0])
for n,w in res:
    a=n.split("_")[0]; arm[a][0]+=w; arm[a][1]+=1
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s\n")
for a in ("POS","FREE","CTRL"):
    if a in arm: print(f"{a:5s} beat champion: {arm[a][0]:4d}/{arm[a][1]:4d} = {arm[a][0]/arm[a][1]:.2%}")
win=[n for n,w in res if w]
json.dump({"winners":win,"seqs":{n:seqs[n] for n in win},"parents":top10},
          open("checkpoints/stage_f.json","w"), indent=1)
print(f"\ntotal winners: {len(win)}")
for n in win[:25]:
    s=seqs[n]; print(f"  {n}  TATA={tata(s)} GC={gc(s):.3f}")
