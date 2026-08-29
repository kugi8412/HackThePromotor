"""40 calls: 20 random pairs from the natural panel, each asked both ways.

Purpose is to measure the b-slot win rate in *this* regime (diverged natural
promoters), because our existing 4% figure comes from near-wild variants and the
Judge is documented not to transfer between regimes.
"""
import random, json
from judge_scan import load_panel, CachedJudge, tata, gc

panel = load_panel()
rng = random.Random(20260829)
pairs = []
seen = set()
while len(pairs) < 20:
    i, j = rng.sample(range(len(panel)), 2)
    if (i, j) in seen or (j, i) in seen: continue
    seen.add((i, j)); pairs.append((i, j))

J = CachedJudge()
rows = []
for i, j in pairs:
    (na, sa), (nb, sb) = panel[i], panel[j]
    fwd = J(na, sa, nb, sb)      # na in slot a
    rev = J(nb, sb, na, sa)      # nb in slot a
    # fwd True  = na won from slot a ; rev True = nb won from slot a
    if fwd and not rev:   verdict = "A_real"     # na won both slots
    elif rev and not fwd: verdict = "B_real"     # nb won both slots
    elif fwd and rev:     verdict = "tie_first"  # each won when placed first
    else:                 verdict = "tie_second" # each LOST when placed first
    rows.append(dict(a=na, b=nb, fwd=fwd, rev=rev, verdict=verdict))
    print(f"{na} vs {nb}: fwd={int(fwd)} rev={int(rev)}  {verdict}")

J.save()
from collections import Counter
c = Counter(r["verdict"] for r in rows)
print("\n--- verdicts over 20 pairs ---")
for k in ("A_real", "B_real", "tie_first", "tie_second"):
    print(f"  {k:12s} {c.get(k,0):2d}")
real = c.get("A_real", 0) + c.get("B_real", 0)
print(f"\nreal preferences: {real}/20 = {real/20:.0%}")
# b-slot win rate across all 40 individual queries
bwins = sum(1 for r in rows if not r["fwd"]) + sum(1 for r in rows if not r["rev"])
print(f"slot-b wins: {bwins}/40 = {bwins/40:.1%}   (near-wild reference was 1/24 = 4.2%)")
print(f"paid calls: {J.paid}  cache hits: {J.hits}")
json.dump(rows, open("checkpoints/calibration_natural.json", "w"), indent=1)
