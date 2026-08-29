#!/usr/bin/env python3
import sys, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
hits_f, meta_f, outp = sys.argv[1:4]
meta = pd.read_csv(meta_f, sep="\t", dtype=str)
hdr, rows = None, []
for l in open(hits_f):
    l = l.rstrip("\n")
    if not l.strip(): continue
    if l.startswith("#"):
        hdr = l.lstrip("#").split("\t")      # format=3 header starts with '#'
        continue
    rows.append(l.split("\t"))
if hdr is None: hdr = rows.pop(0)
hits = pd.DataFrame(rows, columns=[c.strip() for c in hdr])
qcol = next(c for c in hits.columns if c.lower().startswith("query"))
rcol = next(c for c in hits.columns if c.lower().startswith("ref"))
hits["transcript_id"] = hits[qcol].str.split().str[0].str.replace(r"\.\d+$", "", regex=True)
hits["chimp_best_hit"] = hits[rcol].str.split().str[0]
for c in ["ANI","WKID","KID","Complt"]:
    if c in hits.columns: hits[c] = pd.to_numeric(hits[c].str.rstrip("%"), errors="coerce")
if "Complt" not in hits.columns and {"KID","WKID"} <= set(hits.columns):
    hits["Complt"] = (100 * hits.KID / hits.WKID).clip(upper=100)   # BBTools 38.90 format=3 lacks Complt
keep = ["transcript_id","chimp_best_hit"] + [c for c in ["ANI","Complt","WKID"] if c in hits.columns]
df = meta.merge(hits[keep], on="transcript_id", how="left")
df["has_hit"] = df.ANI.notna()
def binf(r):
    if not r.has_hit: return "no_hit"
    if r.ANI >= 98 and r.Complt >= 90: return "highly_conserved"
    if r.ANI >= 95: return "conserved"
    return "divergent"
df["seq_bin"] = df.apply(binf, axis=1)
def coarse(t):
    if t == "protein_coding": return "protein_coding"
    if t == "lncRNA": return "lncRNA"
    if "pseudogene" in t: return "pseudogene"
    if t == "retained_intron": return "retained_intron"
    return "other"
df["type"] = df.transcript_type.map(coarse)
df["stratum"] = df.type.where(df.mane_select != "1", "MANE_Select")
order = ["highly_conserved","conserved","divergent","no_hit"]
cnt = df.groupby(["stratum","seq_bin"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
pct = cnt.div(cnt.sum(axis=1), axis=0).mul(100).round(1)
cnt.to_csv(outp + "_counts.tsv", sep="\t"); pct.to_csv(outp + "_pct.tsv", sep="\t")
df.to_csv(outp + "_per_transcript.tsv", sep="\t", index=False)
print(cnt, "\n", pct)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
pct.plot(kind="barh", stacked=True, ax=ax[0]); ax[0].set_xlabel("% of human transcripts")
sns.histplot(df[df.has_hit], x="ANI", hue="stratum", element="step", bins=60, ax=ax[1])
ax[1].set_xlabel("estimated identity to best chimp transcript (%)")
plt.tight_layout(); plt.savefig(outp + "_stage0.png", dpi=150)
