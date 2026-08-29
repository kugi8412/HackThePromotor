"""Stage D: screen everything against our champion TrAtP1_000297.

Three arms:
  1. the 994 stage-A rejects vs champion -- predicted NULL by the threshold model
     (they already score below wild+t, and the champion scores above it), so this
     is a falsification test, not a search.
  2. the organisers' 100 reference promoters (19 species) vs champion AND vs wild
     -- a genuinely untapped pool, the arm that could actually find something.
  3. the 8 designed constructs vs champion.
"""
import csv, json, threading, time
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge, tata, gc

panel = dict(load_panel())
survivors = set(json.load(open("checkpoints/stage_a_winners.json")))
CHAMP = "TrAtP1_000297"
champ = panel[CHAMP]
wild = open("wild_type.txt").read().strip()
rejects = [n for n in panel if n not in survivors]

ref = {}
for r in csv.DictReader(open("Promotory.csv"), delimiter=";"):
    ref["REF_" + r["nazwa"][:30]] = r["sekwencja"]

p = r"C:\Users\leszk\Downloads\Ta_pks1_promoter_panel_800bp.fasta"
des, nm, buf = {}, None, []
for line in open(p):
    line = line.strip()
    if line.startswith(">"):
        if nm: des[nm] = "".join(buf)
        nm, buf = line[1:].split()[0], []
    else: buf.append(line)
des[nm] = "".join(buf)
des = {k[:36]: (wild[:4]+v if len(v)==796 else v) for k,v in des.items() if not k.startswith("P0")}

jobs = []
for n in rejects:            jobs.append(("rejects_vs_champ", CHAMP, champ, n, panel[n]))
for n,s in ref.items():      jobs.append(("ref_vs_champ",     CHAMP, champ, n, s))
for n,s in ref.items():      jobs.append(("ref_vs_wild",      "WILD", wild,  n, s))
for n,s in des.items():      jobs.append(("designed_vs_champ",CHAMP, champ, n, s))
print(f"{len(jobs)} queries queued", flush=True)

J = CachedJudge(); lock = threading.Lock(); done=[0]; t0=time.time()
def run(job):
    arm, an, aseq, bn, bseq = job
    k = f"{an}|{bn}"
    with lock: hit = J.cache.get(k)
    if hit is None:
        v = J.h.sedzia(aseq, bseq, an, bn)
        with lock: J.cache[k]=v; J.paid+=1
    else:
        v = hit
        with lock: J.hits+=1
    with lock:
        done[0]+=1
        if done[0]%300==0:
            print(f"  {done[0]}/{len(jobs)} paid={J.paid} {done[0]/(time.time()-t0):.0f}/s", flush=True)
            J.save()
    return arm, bn, (not v)      # True = challenger beat the anchor from slot b

with ThreadPoolExecutor(max_workers=12) as ex:
    res = list(ex.map(run, jobs))
J.save()

from collections import defaultdict
byarm = defaultdict(list)
for arm, n, won in res: byarm[arm].append((n, won))
print(f"\npaid {J.paid}  cache hits {J.hits}  {time.time()-t0:.0f}s\n")
out={}
for arm in ("rejects_vs_champ","ref_vs_champ","ref_vs_wild","designed_vs_champ"):
    w=[n for n,x in byarm[arm] if x]; tot=len(byarm[arm])
    out[arm]=w
    print(f"{arm:20s} {len(w):3d}/{tot:4d} beat the anchor")
    if w and arm!="rejects_vs_champ": print("     ", ", ".join(w[:20]))
json.dump(out, open("checkpoints/stage_d.json","w"), indent=1)
if not out["rejects_vs_champ"]:
    print("\nthreshold model HELD: 0/994 rejects beat the champion, exactly as predicted.")
else:
    print(f"\nthreshold model FALSIFIED: {len(out['rejects_vs_champ'])} rejects beat the champion.")
