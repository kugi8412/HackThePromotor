"""Copeland + Bradley-Terry/Elo over the stage-B round-robin, plus a transitivity check."""
import json, math, itertools
import numpy as np
from judge_scan import load_panel, CachedJudge, tata, gc

names = sorted(json.load(open("checkpoints/stage_a_winners.json")))
idx = {n: i for i, n in enumerate(names)}
N = len(names)
cache = CachedJudge().cache
panel = dict(load_panel())

D = np.zeros((N, N))          # D[i,j] = 1 if i beat j from BOTH slots
T = np.zeros((N, N))          # T[i,j] = 1 if the pair tied
missing = 0
for a, b in itertools.combinations(names, 2):
    i, j = idx[a], idx[b]
    fwd, rev = cache.get(f"{a}|{b}"), cache.get(f"{b}|{a}")
    if fwd is None or rev is None: missing += 1; continue
    if fwd and not rev:   D[i, j] = 1          # a won in both slots
    elif rev and not fwd: D[j, i] = 1
    else:                 T[i, j] = T[j, i] = 1
print(f"pairs: {N*(N-1)//2}  decisive: {int(D.sum())}  ties: {int(T.sum()//2)}  missing: {missing}")
print(f"decisive fraction: {D.sum()/(N*(N-1)/2):.1%}   (calibration predicted ~55%)")

# ---- transitivity: directed 3-cycles among decisive edges ----
cyc = int(np.trace(D @ D @ D) // 3)
tri = N*(N-1)*(N-2)//6
p_edge = D.sum()/(N*(N-1)/2)
exp_cyc = tri * (p_edge**3) * 2 / 8 * 4   # rough null: both orientations of a 3-cycle
print(f"\n3-cycles (rock-paper-scissors triples): {cyc} of {tri} triples = {cyc/tri:.3%}")

# ---- Copeland ----
wins = D.sum(1); losses = D.sum(0); ties = T.sum(1)
copeland = wins - losses

# ---- Bradley-Terry (MM / Zermelo), ties as half-wins, lightly regularised ----
W = wins + 0.5*ties
p = np.ones(N)
for _ in range(3000):
    S = p[:, None] + p[None, :]
    np.fill_diagonal(S, np.inf)
    denom = (1.0/S).sum(1) + 1.0/(p + 1.0)     # + virtual average opponent (prior)
    new = (W + 0.5) / denom                     # +0.5 pseudo-win keeps zero-win seqs finite
    new /= np.exp(np.log(new).mean())
    if np.max(np.abs(np.log(new) - np.log(p))) < 1e-10: p = new; break
    p = new
elo = 400*np.log10(p/np.exp(np.log(p).mean())) + 1500

order = np.argsort(-elo)
rows = []
for r, i in enumerate(order, 1):
    n = names[i]
    rows.append(dict(rank=r, name=n, elo=round(float(elo[i]),1),
                     copeland=int(copeland[i]), wins=int(wins[i]),
                     losses=int(losses[i]), ties=int(ties[i]),
                     tata=tata(panel[n]), gc=round(gc(panel[n]),4)))
json.dump(rows, open("checkpoints/stage_b_ranking.json","w"), indent=1)


# rank agreement between Copeland and BT, computed without scipy
def rankdata(x):
    o = np.argsort(-x); r = np.empty(len(x)); r[o] = np.arange(len(x)); return r
rc, rb = rankdata(copeland), rankdata(elo)
rho = np.corrcoef(rc, rb)[0,1]
print(f"Copeland vs BT rank correlation: rho = {rho:.4f}")

print(f"\n=== TOP 15 by BT/Elo ===")
print(f"{'#':>3} {'name':16s} {'elo':>7s} {'W':>4s} {'L':>4s} {'T':>4s} {'cope':>5s} {'tata':>4s} {'gc':>6s}")
for r in rows[:15]:
    print(f"{r['rank']:3d} {r['name']:16s} {r['elo']:7.1f} {r['wins']:4d} {r['losses']:4d} "
          f"{r['ties']:4d} {r['copeland']:5d} {r['tata']:4d} {r['gc']:6.3f}")
print("\n=== BOTTOM 5 ===")
for r in rows[-5:]:
    print(f"{r['rank']:3d} {r['name']:16s} {r['elo']:7.1f} {r['wins']:4d} {r['losses']:4d} "
          f"{r['ties']:4d} {r['copeland']:5d} {r['tata']:4d} {r['gc']:6.3f}")

# does TATA predict Elo across the 116?
t = np.array([r["tata"] for r in rows]); e = np.array([r["elo"] for r in rows])
g = np.array([r["gc"] for r in rows])
print(f"\ncorr(TATA, Elo) = {np.corrcoef(t,e)[0,1]:+.3f}   corr(GC, Elo) = {np.corrcoef(g,e)[0,1]:+.3f}")
print("TOP5:", [r["name"] for r in rows[:5]])
json.dump([r["name"] for r in rows[:5]], open("checkpoints/stage_b_top5.json","w"))
