"""
Client for the Hack The Promoter API (hyppe.futura.foundation), plus an adapter
that plugs it into the Engine interface in hack_the_promoter_kit.ipynb.

Read the API notebook carefully before changing anything here. The important
mismatches with our prepared kit:

  * There is NO judge budget. /sedzia is rate limited (600/min), not capped.
  * There is NO batching. One pair per call.
  * The Navigator is NOT an encode/decode autoencoder. There is no way to get a
    latent vector for an arbitrary sequence, and no way to decode an arbitrary
    vector. It offers two things instead:
        /nawigator/mapa    per-position structure map
        /nawigator/edycje  "here are N variants made by perturbing latent codes"
    So the climb-in-latent-space plan does not apply. `propose()` below replaces
    it: the model's own edit proposer becomes the candidate generator, and our
    Judge/Bradley-Terry/portfolio machinery sits on top unchanged.
  * encode() here is a LOCAL featurisation (k-mer + positional composition). It
    costs no API calls and exists only to give the surrogate an input vector.
    It is NOT the model's latent space -- do not claim otherwise.
"""
import json, os, time, urllib.error, urllib.request

import numpy as np

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/128.0 Safari/537.36")
BASES = "ACGT"


def _load_env(path=".env"):
    """Minimal .env reader -- python-dotenv's find_dotenv() breaks under exec()."""
    env = {}
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env


class Hyppe:
    """Thin HTTP client. Mirrors `wolaj` from the organisers' notebook, plus a
    self-throttle so we stay under the per-minute limit instead of relying on
    being told off with a 429."""

    def __init__(self, url=None, key=None, env_path=".env", min_interval=0.0):
        env = _load_env(env_path)
        self.url = (url or os.environ.get("HYPPE_URL") or env.get("HYPPE_URL")
                    or "https://hyppe.futura.foundation").rstrip("/")
        self.key = key or os.environ.get("HYPPE_TOKEN") or env.get("HYPPE_TOKEN")
        if not self.key:
            raise RuntimeError("no API key: put HYPPE_TOKEN in .env")
        self.min_interval = min_interval          # seconds between calls
        self._last = 0.0
        self.calls = 0                            # /sedzia verdicts served
        self.requests = 0                         # every HTTP request
        self.by_path = {}

    def call(self, path, data=None, timeout=900, tries=6):
        if self.min_interval:
            gap = time.time() - self._last
            if gap < self.min_interval: time.sleep(self.min_interval - gap)
        for n in range(tries):
            r = urllib.request.Request(self.url + path)
            r.add_header("X-API-Key", self.key)
            r.add_header("User-Agent", UA)          # without this Cloudflare 1010s
            if data is not None:
                r.add_header("Content-Type", "application/json")
                r.data = json.dumps(data).encode()
            try:
                with urllib.request.urlopen(r, timeout=timeout) as o:
                    self.requests += 1
                    self._last = time.time()
                    self.by_path[path] = self.by_path.get(path, 0) + 1
                    return o.status, json.loads(o.read().decode())
            except urllib.error.HTTPError as e:
                body = e.read().decode("utf-8", "replace")
                # /wgraj 429 means "5 minute cooldown" -- retrying cannot help
                # 5xx here are Cloudflare/origin blips (530 = tunnel down), not our request
                retryable = (e.code in (500, 502, 503, 504, 520, 521, 522, 523, 524, 530)
                             or (e.code == 429 and path != "/wgraj"))
                if not retryable or n == tries - 1:
                    return e.code, body
                wait = float(e.headers.get("Retry-After") or 0) or 0.4 * 2 ** n
                time.sleep(min(wait, 8.0))
            except urllib.error.URLError as e:
                if n == tries - 1: return 0, f"network: {e}"
                time.sleep(0.4 * 2 ** n)
        return 0, "exhausted retries"

    # ---------- endpoints ----------
    def me(self):
        c, r = self.call("/me")
        if c != 200: raise RuntimeError(f"/me -> {c}: {str(r)[:200]}")
        return r

    def dziki(self):
        """The wild-type pks1 promoter: 800 bp, aligned to the TSS."""
        c, r = self.call("/dziki")
        if c != 200: raise RuntimeError(f"/dziki -> {c}: {str(r)[:200]}")
        return r

    def sedzia(self, a, b, nazwa_a="a", nazwa_b="b"):
        """-> True if `a` is the stronger sequence. Ties resolve to whichever the
        server names; treat a single verdict as one noisy observation."""
        c, r = self.call("/sedzia", {"a": a, "b": b,
                                     "nazwa_a": nazwa_a, "nazwa_b": nazwa_b})
        if c != 200: raise RuntimeError(f"/sedzia -> {c}: {str(r)[:200]}")
        self.calls += 1
        return r["silniejsza_idx"] == 0

    def mapa(self, seq, od=0, ile=800):
        c, r = self.call("/nawigator/mapa", {"sekwencja": seq, "od": od, "ile": ile})
        if c != 200: raise RuntimeError(f"/nawigator/mapa -> {c}: {str(r)[:200]}")
        return r

    def edycje(self, seq, poziom=2, ile_kodow=8, opcji=8, ziarno=None):
        d = {"sekwencja": seq, "poziom": poziom, "ile_kodow": ile_kodow, "opcji": opcji}
        if ziarno is not None: d["ziarno"] = ziarno
        c, r = self.call("/nawigator/edycje", d)
        if c != 200: raise RuntimeError(f"/nawigator/edycje -> {c}: {str(r)[:200]}")
        return r

    def wgraj(self, fasta):
        return self.call("/wgraj", {"fasta": fasta})

    def ranking(self):
        c, r = self.call("/ranking")
        if c != 200: raise RuntimeError(f"/ranking -> {c}: {str(r)[:200]}")
        return r


# --------------------------------------------------------------------------
# Backend adapter: satisfies the `Engine` backend contract in the kit.
# --------------------------------------------------------------------------
def kmer_features(seqs, k=4, blocks=8):
    """Local, deterministic featurisation: k-mer frequencies within `blocks`
    positional bins. NOT the model's latent space -- purely a surrogate input,
    and it costs zero API calls."""
    seqs = list(seqs)
    idx = {"".join(p): i for i, p in
           enumerate(__import__("itertools").product(BASES, repeat=k))}
    d = len(idx)
    X = np.zeros((len(seqs), blocks * d), dtype=np.float32)
    for j, s in enumerate(seqs):
        L = len(s); step = max(1, L // blocks)
        for b in range(blocks):
            seg = s[b * step:(b + 1) * step] if b < blocks - 1 else s[b * step:]
            for i in range(len(seg) - k + 1):
                t = idx.get(seg[i:i + k])
                if t is not None: X[j, b * d + t] += 1.0
        n = np.linalg.norm(X[j])
        if n: X[j] /= n
    return X


class HyppeBackend:
    """Drop-in backend for Engine.

    judge()      -> /sedzia
    encode()     -> LOCAL k-mer features (no API); only feeds the surrogate
    propose()    -> /nawigator/edycje; replaces climb()+decode()
    NO judge_batch: the API takes one pair per call, and Engine checks hasattr,
    so the method is deliberately absent rather than a raising stub.
    """
    name = "hyppe"

    def __init__(self, client=None, poziom=2, ile_kodow=8, opcji=8):
        self.h = client or Hyppe()
        self.poziom, self.ile_kodow, self.opcji = poziom, ile_kodow, opcji
        self._seed = 1000

    @property
    def calls(self): return self.h.calls

    def judge(self, a, b):
        return self.h.sedzia(a, b)

    def encode(self, seqs):
        return kmer_features(seqs)

    def decode(self, Z):
        raise NotImplementedError(
            "The Navigator has no decoder for arbitrary vectors. Use propose().")

    def propose(self, seeds, n_out, rng=None):
        """Ask the Navigator for variants of each seed. Returns unique strings."""
        out = {}
        for s in seeds:
            if len(out) >= n_out: break
            self._seed += 1
            try:
                e = self.h.edycje(s, self.poziom, self.ile_kodow, self.opcji, self._seed)
            except RuntimeError as ex:
                print("  edycje failed:", str(ex)[:120]); break
            for o in e.get("opcje", []):
                q = o.get("sekwencja")
                if q and len(q) == 800: out[q] = o.get("nr")
        return list(out)[:n_out]


if __name__ == "__main__":
    h = Hyppe()
    me = h.me()
    print("team      :", me.get("druzyna"), "|", me.get("uczestnik"))
    print("limits/min:", me.get("limity_na_minute"))
    print("upload in :", me.get("zgloszenie_mozliwe_za_s"), "s")
    d = h.dziki()
    print("\nwild type :", d["nazwa"], "| gene", d["gen"], "|", d["genom"])
    print("            %d bp | sha256[:12] %s" % (d["dlugosc"], d["sha256_12"]))
    open("wild_type.txt", "w").write(d["sekwencja"])
    print("            saved -> wild_type.txt")
