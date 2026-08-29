"""Stage J: do sequences that all tie with the champion also tie with EACH OTHER?

Under a threshold model tie is not transitive: |A-champ|<t and |B-champ|<t allows
|A-B| up to 2t, which can exceed t. So decisive pairs are predicted to exist. Zero
decisive pairs across the whole round-robin would instead mean true saturation.
"""
import json, itertools, random, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge

panel=dict(load_panel())
rank=json.load(open("checkpoints/stage_b_ranking.json"))
CHAMP=rank[0]["name"]; plateau=[r["name"] for r in rank if r["losses"]==0]
rng=random.Random(11)

# rebuild the generated sequences we need
top5=[r["name"] for r in rank[:5]]; top10=[r["name"] for r in rank[:10]]
def blocks(s,k):
    b=[round(800*i/k) for i in range(k+1)]; return [s[b[i]:b[i+1]] for i in range(k)]
chim={}
for k in (2,3,4):
    parts={p:blocks(panel[p],k) for p in top5}
    for combo in itertools.product(range(5),repeat=k):
        if len(set(combo))==1: continue
        chim[f"C{k}_"+"".join(map(str,combo))]="".join(parts[top5[c]][i] for i,c in enumerate(combo))
verd=json.load(open("checkpoints/stage_e_verdicts.json"))
ties=[n for n,v in verd.items() if v=="tie"]
pick={}
for n in rng.sample(ties,25): pick[n]=chim[n]

B=10; NB=80
frag={p:[panel[p][i*B:(i+1)*B] for i in range(NB)] for p in top10}
pool=[frag[p][i] for p in top10 for i in range(NB)]
r2=random.Random(20260829)
posl,freel=[],[]
for j in range(3000):
    combo=[r2.randrange(10) for _ in range(NB)]
    posl.append((f"POS_{j:04d}","".join(frag[top10[c]][i] for i,c in enumerate(combo))))
for j in range(3000):
    freel.append((f"FREE_{j:04d}","".join(r2.sample(pool,NB))))
for n,s in rng.sample(posl,15): pick[n]=s
for n,s in rng.sample(freel,5):  pick[n]=s

champ=panel[CHAMP]
for i in rng.sample(range(800),10):
    b=rng.choice([c for c in "ACGT" if c!=champ[i]])
    pick[f"MUT_{i:03d}{b}"]=champ[:i]+b+champ[i+1:]
pick[CHAMP[:36]]=champ
for p in rng.sample([x for x in plateau if x!=CHAMP],4): pick[p]=panel[p]

names=list(pick)
print(f"{len(names)} sequences -> {len(names)*(len(names)-1)} ordered queries", flush=True)
J=CachedJudge(); lock=threading.Lock(); done=[0]; t0=time.time()
def run(pr):
    a,b=pr; key=f"{a}|{b}"
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(pick[a],pick[b],a[:36],b[:36])
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%1000==0: print(f"  {done[0]} paid={J.paid}",flush=True); J.save()
    return key,v
pairs=list(itertools.permutations(names,2))
list(ThreadPoolExecutor(max_workers=12).map(run,pairs))
J.save()
c=J.cache; dec=[]; tie=0
for a,b in itertools.combinations(names,2):
    f,r=c[f"{a}|{b}"],c[f"{b}|{a}"]
    if f and not r: dec.append((a,b))
    elif r and not f: dec.append((b,a))
    else: tie+=1
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print(f"pairs {len(dec)+tie}  decisive {len(dec)}  ties {tie}  ({len(dec)/(len(dec)+tie):.1%} decisive)")
if dec:
    from collections import Counter
    w=Counter(a for a,_ in dec); l=Counter(b for _,b in dec)
    print("\ntop winners inside the tie-set:")
    for n,k in w.most_common(12): print(f"  {n:22s} beats {k:3d}  loses {l[n]:3d}")
json.dump({"decisive":dec,"names":names}, open("checkpoints/stage_j.json","w"), indent=1)
