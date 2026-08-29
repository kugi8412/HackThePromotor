"""Stage A: screen all 1110 natural promoters against the wild-type pks1 promoter.

Wild type sits in slot `a` (the advantaged slot) for every query, so a candidate
that wins has won from the disadvantaged slot -- a trustworthy win. Using a fixed
bar rather than a rising champion keeps all 1110 measurements comparable.
"""
import json, gzip, os, threading, time
from concurrent.futures import ThreadPoolExecutor
import hyppe
from judge_scan import load_panel, CachedJudge, tata, gc

panel = load_panel()
wild = open("wild_type.txt").read().strip()
assert len(wild) == 800

J = CachedJudge()
lock = threading.Lock()
done = [0]
t0 = time.time()

def screen(item):
    name, seq = item
    with lock:
        k = f"WILD|{name}"
        hit = J.cache.get(k)
    if hit is not None:
        wild_won = hit
    else:
        wild_won = J.h.sedzia(wild, seq, "WILD", name)   # wild in slot a
        with lock:
            J.cache[k] = wild_won
            J.paid += 1
    with lock:
        done[0] += 1
        if done[0] % 200 == 0:
            el = time.time() - t0
            print(f"  {done[0]}/{len(panel)}  {done[0]/el:.0f} calls/s  paid={J.paid}", flush=True)
    return name, (not wild_won)      # True = candidate beat wild from slot b

with ThreadPoolExecutor(max_workers=12) as ex:
    results = list(ex.map(screen, panel))

J.save()
seqs = dict(panel)
winners = [n for n, w in results if w]
rows = [dict(name=n, beat_wild=w, tata=tata(seqs[n]), gc=round(gc(seqs[n]), 4))
        for n, w in results]
json.dump(rows, open("checkpoints/stage_a_natural.json", "w"), indent=1)

print(f"\n=== STAGE A ===")
print(f"paid calls: {J.paid}   cache hits: {J.hits}   elapsed {time.time()-t0:.0f}s")
print(f"beat wild type from slot b: {len(winners)}/{len(panel)} = {len(winners)/len(panel):.1%}")

# TATA confound check
import statistics
w_t = [r["tata"] for r in rows if r["beat_wild"]]
l_t = [r["tata"] for r in rows if not r["beat_wild"]]
w_g = [r["gc"] for r in rows if r["beat_wild"]]
l_g = [r["gc"] for r in rows if not r["beat_wild"]]
print(f"\nTATA per 800bp:  winners {statistics.mean(w_t):.3f} (n={len(w_t)})   "
      f"losers {statistics.mean(l_t):.3f} (n={len(l_t)})")
print(f"  TATA-free:     winners {sum(1 for t in w_t if t==0)/len(w_t):.1%}   "
      f"losers {sum(1 for t in l_t if t==0)/len(l_t):.1%}")
print(f"GC:              winners {statistics.mean(w_g):.4f}   losers {statistics.mean(l_g):.4f}")
json.dump(winners, open("checkpoints/stage_a_winners.json", "w"), indent=1)
