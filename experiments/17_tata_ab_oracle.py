"""A/B the Oracle on TATA status, using the per-upload reported positions.

best-of only reports "better than your own best", which both files will fail, so the
readout is `pozycja_top100` from each upload response. /ranking is snapshotted around
each upload: if the field moves between them the comparison is void.
"""
import hyppe, json, time
h=hyppe.Hyppe()
def field():
    r=h.ranking()
    return {e["druzyna"]:(e["punkty_top10"],e["punkty_top100"],e["wgranie_o"]) for e in r["ranking"]}
def up(fn):
    while True:
        s=h.me()["zgloszenie_mozliwe_za_s"]
        if s<=0: break
        print(f"  cooldown {s}s ...", flush=True); time.sleep(min(s+2,320))
    c,r=h.wgraj(open(fn).read())
    return c,r
res={}
f0=field()
for tag,fn in (("L_tataless","submission_groupL_tataless.fasta"),
               ("T_tata","submission_groupT_tata.fasta")):
    before=field()
    c,r=up(fn)
    after=field()
    if c!=200:
        print(tag,"HTTP",c,str(r)[:200]); continue
    res[tag]=dict(top10=r["pozycja_top10"], top100=r["pozycja_top100"],
                  punkty=r["punkty_razem"], ocen=r["ocenionych"])
    others_moved = {k:v for k,v in after.items() if k!="druzyna_02" and before.get(k)!=v}
    res[tag]["field_moved"]=others_moved
    print(f"{tag}: evaluated={r['ocenionych']} top10_pos={r['pozycja_top10']} "
          f"top100_pos={r['pozycja_top100']} punkty={r['punkty_razem']}"
          f"{'  [FIELD MOVED: '+str(list(others_moved))+']' if others_moved else ''}", flush=True)
f1=field()
moved={k:(f0.get(k),f1.get(k)) for k in set(f0)|set(f1) if k!="druzyna_02" and f0.get(k)!=f1.get(k)}
print("\n=== RESULT ===")
for k,v in res.items(): print(" ",k,v)
print("\nfield changes across the whole experiment (excluding us):", moved if moved else "NONE - comparison valid")
print("our banked best still:", f1.get("druzyna_02"))
json.dump({"res":res,"field_before":f0,"field_after":f1,"moved":moved}, open("checkpoints/tata_ab.json","w"))
