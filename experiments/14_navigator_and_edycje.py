"""Stage L: (a) Navigator-safety top-20 single mutants, (b) edycje variants screened.

(a) ranks the 2,396 champion-tying mutants by SAFETY, not superiority -- none is
    better, so the criterion is Navigator endorsement first, then low wagaP and
    rekon=0, avoiding ->G in the high-leverage 701-800 tail.
(b) pulls latent-code variants off the champion and screens them against it.
"""
import json, time, threading, textwrap
from concurrent.futures import ThreadPoolExecutor
import hyppe
from judge_scan import load_panel, CachedJudge, tata, gc

panel=dict(load_panel()); rank=json.load(open("checkpoints/stage_b_ranking.json"))
CH=rank[0]["name"]; champ=panel[CH]
mp=json.load(open("checkpoints/champ_mapa.json"))["pozycje"]
verd=json.load(open("checkpoints/stage_k_mutant_verdicts.json"))
byidx={p["poz"]-1:p for p in mp}

# ---------- (a) safety ranking ----------
cand=[]
for n,v in verd.items():
    if v!="tie": continue
    i=int(n[3:6]); b=n[6]; p=byidx[i]
    endorsed = (p["zmien_na"]==b)
    risky = (b=="G" and i>=700)
    cand.append(dict(name=n,pos=i,base=b,wagaP=p["wagaP"],rekon=p["rekon"],
                     endorsed=endorsed,risky=risky))
cand.sort(key=lambda c:(not c["endorsed"], c["risky"], c["wagaP"], c["rekon"]))
top20=cand[:20]
print("=== (a) Navigator-safety top 20 single mutants ===")
print(f"{'#':>3} {'mutation':12s} {'TSS':>6s} {'wagaP':>6s} {'rekon':>5s} {'endorsed':>8s}")
for j,c in enumerate(top20,1):
    print(f"{j:3d} {c['base']+'@'+str(c['pos']+1):12s} {c['pos']+1-800:6d} {c['wagaP']:6.3f} {c['rekon']:5d} {str(c['endorsed']):>8s}")
with open("top20_champion_mutants.fasta","w") as f:
    for j,c in enumerate(top20,1):
        s=champ[:c["pos"]]+c["base"]+champ[c["pos"]+1:]
        f.write(f">{CH}_mut{j:02d}|pos={c['pos']+1}|TSS{c['pos']+1-800:+d}|{champ[c['pos']]}to{c['base']}"
                f"|wagaP={c['wagaP']}|navigator_endorsed={int(c['endorsed'])}|judge=tie\n")
        f.write("\n".join(textwrap.wrap(s,80))+"\n")
print("wrote top20_champion_mutants.fasta")

# ---------- (b) edycje pool ----------
h=hyppe.Hyppe()
def edy(poz, seed):
    for a in range(6):
        try: return h.edycje(champ, poziom=poz, ile_kodow=8, opcji=8, ziarno=seed)
        except Exception: time.sleep(3)
    return None
pool={}
lock=threading.Lock()
jobs=[(p,s) for p in (1,2,3) for s in range(1,61)]
def gen(j):
    p,s=j; e=edy(p,s)
    if not e: return
    with lock:
        for o in e["opcje"]:
            pool[f"E{p}_{s:02d}_{o['nr']}"]=o["sekwencja"]
list(ThreadPoolExecutor(max_workers=8).map(gen, jobs))
pool={k:v for k,v in pool.items() if len(v)==800 and set(v)<=set("ACGT")}
uniq={}
for k,v in pool.items(): uniq.setdefault(v,k)
pool={v:k for k,v in uniq.items()}
print(f"\n=== (b) edycje pool: {len(pool)} unique variants ===")
d=[sum(1 for a,b in zip(s,champ) if a!=b) for s in pool.values()]
t=[tata(s) for s in pool.values()]
print(f"distance to champion: min {min(d)} median {sorted(d)[len(d)//2]} max {max(d)}")
print(f"TATA boxes: mean {sum(t)/len(t):.2f} (champion has {tata(champ)});  {sum(1 for x in t if x>0)}/{len(t)} carry >=1")

J=CachedJudge(); done=[0]; t0=time.time()
def run(it):
    n,s=it; key=f"{CH[:36]}|{n}"
    with lock: hit=J.cache.get(key)
    if hit is None:
        v=J.h.sedzia(champ,s,CH[:36],n[:36])
        with lock: J.cache[key]=v; J.paid+=1
    else:
        v=hit
        with lock: J.hits+=1
    key2=f"{n}|{CH[:36]}"
    with lock: hit2=J.cache.get(key2)
    if hit2 is None:
        v2=J.h.sedzia(s,champ,n[:36],CH[:36])
        with lock: J.cache[key2]=v2; J.paid+=1
    else:
        v2=hit2
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%300==0: print(f"  {done[0]}/{len(pool)} paid={J.paid}",flush=True); J.save()
    return n,("BEATS" if (not v and v2) else "loses" if (v and not v2) else "tie")
res=list(ThreadPoolExecutor(max_workers=12).map(run, pool.items()))
J.save()
from collections import Counter
print(f"\npaid {J.paid} cache hits {J.hits} {time.time()-t0:.0f}s")
print("verdict vs champion:", dict(Counter(v for _,v in res)))
win=[n for n,v in res if v=="BEATS"]
print("winners:", win[:20] if win else "none")
json.dump({"verdicts":dict(res),"seqs":{n:pool[n] for n,_ in res}}, open("checkpoints/stage_l.json","w"))
