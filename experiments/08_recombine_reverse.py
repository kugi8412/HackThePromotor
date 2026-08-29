"""Stage E reverse: chimera in slot `a`, champion in slot `b`.

Combined with the forward screen this gives a three-way verdict per chimera.
Champion winning from the disadvantaged slot means the chimera is strictly worse;
chimera winning from its advantaged slot, having failed to win from slot b, means
the pair is a tie.
"""
import json, itertools, threading, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge

panel=dict(load_panel()); top5=json.load(open("checkpoints/stage_b_top5.json"))
CHAMP=top5[0]; champ=panel[CHAMP]
def blocks(s,k):
    b=[round(800*i/k) for i in range(k+1)]; return [s[b[i]:b[i+1]] for i in range(k)]
chim={}
for k in (2,3,4):
    parts={p:blocks(panel[p],k) for p in top5}
    for combo in itertools.product(range(5),repeat=k):
        if len(set(combo))==1: continue
        s="".join(parts[top5[c]][i] for i,c in enumerate(combo))
        chim.setdefault(s,(k,combo))
items=[(f"C{k}_"+"".join(map(str,c)), s, k, c) for s,(k,c) in chim.items()]
print(f"{len(items)} reverse queries", flush=True)

J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(it):
    n,s,k,c = it
    key=f"{n}|{CHAMP}"
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(s, champ, n[:36], CHAMP[:36])
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%200==0: print(f"  {done[0]}/{len(items)} paid={J.paid}",flush=True); J.save()
    return n,k,c,v      # v True = chimera won from slot a

res=list(ThreadPoolExecutor(max_workers=12).map(run, items))
J.save()
cache=J.cache
byk=defaultdict(lambda: defaultdict(int)); verdicts={}
for n,k,c,rev in res:
    fwd = cache[f"{CHAMP}|{n}"]          # True = champion won from slot a
    if not fwd and rev:  v="chimera_wins"
    elif fwd and not rev: v="LOSES_to_champion"
    else:                 v="tie"
    verdicts[n]=v; byk[k][v]+=1
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print(f"\n{'k':>2} {'n':>4} {'tie':>6} {'loses':>6} {'wins':>5}  loss rate")
for k in sorted(byk):
    tot=sum(byk[k].values()); l=byk[k]["LOSES_to_champion"]
    print(f"{k:2d} {tot:4d} {byk[k]['tie']:6d} {l:6d} {byk[k]['chimera_wins']:5d}  {l/tot:6.1%}")
tot=len(res); L=sum(1 for v in verdicts.values() if v=="LOSES_to_champion")
print(f"\noverall: {L}/{tot} = {L/tot:.1%} strictly worse than the champion")
json.dump(verdicts, open("checkpoints/stage_e_verdicts.json","w"), indent=1)
