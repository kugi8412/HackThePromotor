"""Stage E: recombine the top-5 Elo naturals and screen the chimeras against the champion.

Split each parent into k equal blocks (k = 2,3,4) and take every combination of
parental origins per block. Parents and cross-k duplicates (a 4-way AABB is the
same sequence as a 2-way AB) are removed. Champion sits in slot `a`, so any
chimera that wins has won from the disadvantaged slot.
"""
import json, itertools, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge, tata, gc

panel = dict(load_panel())
top5 = json.load(open("checkpoints/stage_b_top5.json"))
CHAMP = top5[0]; champ = panel[CHAMP]
print("parents:", top5)

def blocks(seq, k):
    b = [round(800*i/k) for i in range(k+1)]
    return [seq[b[i]:b[i+1]] for i in range(k)]

chim = {}
for k in (2,3,4):
    parts = {p: blocks(panel[p], k) for p in top5}
    for combo in itertools.product(range(5), repeat=k):
        if len(set(combo)) == 1: continue                 # that's just a parent
        s = "".join(parts[top5[c]][i] for i,c in enumerate(combo))
        assert len(s) == 800
        chim.setdefault(s, f"C{k}_" + "".join(str(c) for c in combo))
seqs = {v:k for k,v in chim.items()}
# parents as controls: all tie with the champion in the round-robin, so should tie here
for p in top5[1:]: seqs["CTRL_"+p[-6:]] = panel[p]
print(f"unique chimeras: {len(chim)}  (+{len(top5)-1} parent controls) -> {len(seqs)} calls")

J = CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(item):
    n, s = item
    k = f"{CHAMP}|{n}"
    with lock: hit = J.cache.get(k)
    if hit is None:
        v = J.h.sedzia(champ, s, CHAMP[:36], n[:36])
        with lock: J.cache[k]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%200==0: print(f"  {done[0]}/{len(seqs)} paid={J.paid}", flush=True); J.save()
    return n, (not v)          # True = chimera beat the champion from slot b

with ThreadPoolExecutor(max_workers=12) as ex:
    res = list(ex.map(run, seqs.items()))
J.save()
win = [n for n,w in res if w]
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print(f"chimeras beating the champion: {len(win)}/{len(seqs)}")
ctrl = [n for n,w in res if n.startswith('CTRL_')]
print(f"parent controls that beat champion: {sum(1 for n,w in res if n.startswith('CTRL_') and w)}/{len(ctrl)}  (expect 0)")

# does the LAST block's donor decide the outcome?
from collections import defaultdict
tail = defaultdict(lambda:[0,0])
for n,w in res:
    if not n.startswith("C"): continue
    d = n.split("_")[1][-1]                # donor of the final block
    tail[d][0]+=w; tail[d][1]+=1
print("\nwin rate by donor of the FINAL block (index into top5):")
for d in sorted(tail): print(f"  parent {d} ({top5[int(d)]}): {tail[d][0]}/{tail[d][1]} = {tail[d][0]/tail[d][1]:.1%}")
head = defaultdict(lambda:[0,0])
for n,w in res:
    if not n.startswith("C"): continue
    d = n.split("_")[1][0]
    head[d][0]+=w; head[d][1]+=1
print("win rate by donor of the FIRST block:")
for d in sorted(head): print(f"  parent {d} ({top5[int(d)]}): {head[d][0]}/{head[d][1]} = {head[d][0]/head[d][1]:.1%}")
json.dump({"winners":win,"seqs":{n:seqs[n] for n in win}}, open("checkpoints/stage_e.json","w"), indent=1)
