"""Stage C: the top 5 BT-ranked natural promoters against the designed panel.

Both orders, so each of the 40 matchups yields a three-way verdict rather than a
positional artefact. P6/P8 are padded to 800 bp at the 5' distal end, as in the
earlier panel test.
"""
import json, itertools
from judge_scan import load_panel, CachedJudge, tata, gc

top5 = json.load(open("checkpoints/stage_b_top5.json"))
panel = dict(load_panel())
wild = open("wild_type.txt").read().strip()
pad = wild[:4]

p = r"C:\Users\leszk\Downloads\Ta_pks1_promoter_panel_800bp.fasta"
des, name, buf = {}, None, []
for line in open(p):
    line = line.strip()
    if line.startswith(">"):
        if name: des[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else: buf.append(line)
des[name] = "".join(buf)
des = {k: (pad + v if len(v) == 796 else v) for k, v in des.items()}
designed = [k for k in des if not k.startswith("P0")]      # exclude the pks1 baselines

J = CachedJudge()
rows = []
print(f"{'natural':16s} {'designed':40s} {'verdict'}")
for nat in top5:
    for dz in designed:
        a, b = nat, dz[:36]
        fwd = J(a, panel[nat], b, des[dz])     # natural in slot a
        rev = J(b, des[dz], a, panel[nat])     # designed in slot a
        if fwd and not rev:   v = "natural_wins"
        elif rev and not fwd: v = "designed_wins"
        else:                 v = "tie"
        rows.append(dict(natural=nat, designed=dz, fwd=fwd, rev=rev, verdict=v))
        print(f"{nat:16s} {dz:40s} {v}")
J.save()
json.dump(rows, open("checkpoints/stage_c_top5_vs_designed.json","w"), indent=1)

from collections import Counter
c = Counter(r["verdict"] for r in rows)
print(f"\n=== {len(rows)} matchups ===")
for k in ("natural_wins","designed_wins","tie"): print(f"  {k:15s} {c.get(k,0)}")
print("\nby designed construct:")
for dz in designed:
    sub = [r for r in rows if r["designed"] == dz]
    cc = Counter(r["verdict"] for r in sub)
    print(f"  {dz:40s} wins {cc.get('designed_wins',0)}/5  ties {cc.get('tie',0)}  losses {cc.get('natural_wins',0)}")
print("\nby natural promoter:")
for nat in top5:
    sub = [r for r in rows if r["natural"] == nat]
    cc = Counter(r["verdict"] for r in sub)
    print(f"  {nat:16s} wins {cc.get('natural_wins',0)}/{len(sub)}  ties {cc.get('tie',0)}  losses {cc.get('designed_wins',0)}")
print(f"\npaid {J.paid}  cache hits {J.hits}")
