"""Screening the T. atroviride natural-promoter panel against the hyppe Judge.

Protocol notes live in memory; the short version is that /sedzia is deterministic
and near-totally biased toward slot `a`, so a single verdict mostly reports which
sequence was listed first. Everything here is built on the asymmetric screen: put
the incumbent in slot `a`, and treat only a win from slot `b` as real.
"""
import gzip, json, os, random, re, sys, time
import hyppe

FASTA = r"C:\Users\leszk\Downloads\T_atroviride_promoters.fasta"
CACHE = "checkpoints/natural_scan_cache.json.gz"


def load_panel(path=FASTA, keep=800):
    """Records are 1000 bp and ATG-anchored; the last `keep` bp is the ATG-proximal
    end, which is the closest we get to the organisers' TSS frame."""
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name: out.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name: out.append((name, "".join(buf)))
    return [(n.split("|")[0], s[-keep:]) for n, s in out]


def tata(s):
    return len(re.findall(r"(?=TATA[AT]A[AT])", s))


def gc(s):
    return (s.count("G") + s.count("C")) / len(s)


class CachedJudge:
    """Wraps /sedzia with an on-disk cache. The Judge is deterministic, so a repeat
    query buys nothing -- every cache hit is a call we do not have to pay for."""

    def __init__(self, client=None, path=CACHE):
        self.h = client or hyppe.Hyppe()
        self.path = path
        self.cache = {}
        if os.path.exists(path):
            self.cache = json.load(gzip.open(path))
        self.paid = 0
        self.hits = 0

    def __call__(self, ida, seqa, idb, seqb):
        """-> True if the sequence in slot `a` won."""
        k = f"{ida}|{idb}"
        if k in self.cache:
            self.hits += 1
            return self.cache[k]
        v = self.h.sedzia(seqa, seqb, ida, idb)
        self.cache[k] = v
        self.paid += 1
        if self.paid % 250 == 0: self.save()
        return v

    def save(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with gzip.open(self.path, "wt") as f:
            json.dump(self.cache, f)
