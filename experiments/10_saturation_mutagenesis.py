"""Stage G: mutational scanning of the top-5 naturals against the champion.

Arm 1 -- exhaustive single-nucleotide saturation scan (800 positions x 3 alternative
         bases x 5 parents). Predicted null: 1 bp is far below the measured ~25 bp
         resolution floor, so nothing can cross the Judge's threshold.
Arm 2 -- window-swap scan, pitched ABOVE that floor: one 50 or 100 bp window of the
         champion replaced by the corresponding window of another top-5 parent.
"""
import json, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge

panel=dict(load_panel())
rank=json.load(open("checkpoints/stage_b_ranking.json"))
top5=[r["name"] for r in rank[:5]]
CHAMP=top5[0]; champ=panel[CHAMP]

seqs={}
# arm 1: single-nucleotide saturation
for pi,p in enumerate(top5):
    s=panel[p]
    for i in range(800):
        for base in "ACGT":
            if base==s[i]: continue
            seqs[f"M{pi}_{i:03d}{base}"]=s[:i]+base+s[i+1:]
n1=len(seqs)
# arm 2: window swaps into the champion, above the resolution floor
for W in (50,100):
    for w in range(800//W):
        for di,d in enumerate(top5[1:],1):
            s=champ[:w*W]+panel[d][w*W:(w+1)*W]+champ[(w+1)*W:]
            seqs[f"W{W}_{w:02d}_{di}"]=s
n2=len(seqs)-n1
assert all(len(v)==800 for v in seqs.values())
print(f"arm1 single-nt: {n1}   arm2 window-swap: {n2}   total {len(seqs)}", flush=True)

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
        if done[0]%2000==0: print(f"  {done[0]}/{len(seqs)} paid={J.paid} {done[0]/(time.time()-t0):.0f}/s",flush=True); J.save()
    return n,(not v)
res=list(ThreadPoolExecutor(max_workers=12).map(run, seqs.items()))
J.save()
a1=[(n,w) for n,w in res if n.startswith("M")]
a2=[(n,w) for n,w in res if n.startswith("W")]
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print(f"arm1 single-nt  beat champion: {sum(w for _,w in a1)}/{len(a1)}")
print(f"arm2 window-swap beat champion: {sum(w for _,w in a2)}/{len(a2)}")
win=[n for n,w in res if w]
if win:
    print("\nWINNERS:")
    for n in win[:40]: print("  ",n)
json.dump({"winners":win,"seqs":{n:seqs[n] for n in win}}, open("checkpoints/stage_g.json","w"), indent=1)
