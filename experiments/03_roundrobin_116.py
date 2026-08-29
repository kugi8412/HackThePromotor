"""Stage B: exhaustive both-orders round-robin over the 116 stage-A survivors.

Every unordered pair is asked twice, once in each slot order. That yields a
three-way verdict per pair -- decisive win, or tie -- instead of the coin flip a
single randomised query gives against a positionally biased judge. Bradley-Terry
fitted on these verdicts is legitimate; fitted on single randomised queries it was
not, which is what invalidated our earlier tournament.
"""
import json, threading, time, itertools
from concurrent.futures import ThreadPoolExecutor
from judge_scan import load_panel, CachedJudge

names = sorted(json.load(open("checkpoints/stage_a_winners.json")))
panel = dict(load_panel())
seqs = {n: panel[n] for n in names}
N = len(names)
pairs = [(a, b) for a, b in itertools.permutations(names, 2)]   # ordered
print(f"{N} sequences -> {len(pairs)} ordered queries", flush=True)

J = CachedJudge()
lock = threading.Lock()
done = [0]; t0 = time.time()

def ask(pr):
    a, b = pr
    k = f"{a}|{b}"
    with lock:
        hit = J.cache.get(k)
    if hit is None:
        v = J.h.sedzia(seqs[a], seqs[b], a, b)      # True = a (slot a) won
        with lock:
            J.cache[k] = v; J.paid += 1
    else:
        v = hit
        with lock: J.hits += 1
    with lock:
        done[0] += 1
        if done[0] % 2000 == 0:
            el = time.time() - t0
            print(f"  {done[0]}/{len(pairs)}  {done[0]/el:.0f}/s  paid={J.paid}", flush=True)
            J.save()
    return k, v

with ThreadPoolExecutor(max_workers=12) as ex:
    for _ in ex.map(ask, pairs): pass
J.save()
print(f"done: paid {J.paid}, cache hits {J.hits}, {time.time()-t0:.0f}s", flush=True)
