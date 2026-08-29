"""Both-orders test of the designed panel against wild-type pks1.

Both orders (2 calls/sequence) gives a three-way verdict: a sequence that wins
from both slots is a real preference; one that wins only when listed first is a
tie. P0 is a positive control -- 1 SNP from wild, far below the ~25bp resolution
floor, so it must come back a tie.
"""
import json
import hyppe
from judge_scan import CachedJudge, tata, gc

p = r"C:\Users\leszk\Downloads\Ta_pks1_promoter_panel_800bp.fasta"
recs, name, buf = [], None, []
for line in open(p):
    line = line.strip()
    if line.startswith(">"):
        if name: recs.append((name, "".join(buf)))
        name, buf = line[1:].split()[0], []
    else: buf.append(line)
recs.append((name, "".join(buf)))

wild = open("wild_type.txt").read().strip()
J = CachedJudge()
rows = []
for n, s in recs:
    try:
        fwd = J("WILD", wild, n, s)     # wild in slot a
        rev = J(n, s, "WILD", wild)     # candidate in slot a
    except RuntimeError as e:
        print(f"{n:40s} ERROR: {str(e)[:120]}")
        rows.append(dict(name=n, verdict="error", err=str(e)[:200])); continue
    if not fwd and rev:   v = "BEATS_WILD"      # won from both slots
    elif fwd and not rev: v = "loses_to_wild"   # wild won from both slots
    elif fwd and rev:     v = "tie"             # each wins when listed first
    else:                 v = "tie_inverted"    # each loses when listed first
    rows.append(dict(name=n, fwd=fwd, rev=rev, verdict=v, len=len(s),
                     tata=tata(s), gc=round(gc(s), 4)))
    print(f"{n:40s} len={len(s)} fwd={int(fwd)} rev={int(rev)}  {v}")

J.save()
json.dump(rows, open("checkpoints/panel_vs_wild.json", "w"), indent=1)
print(f"\npaid {J.paid}  cache hits {J.hits}")
print("\n=== BEATS WILD TYPE (won from both slots) ===")
for r in rows:
    if r.get("verdict") == "BEATS_WILD": print("  ", r["name"])
