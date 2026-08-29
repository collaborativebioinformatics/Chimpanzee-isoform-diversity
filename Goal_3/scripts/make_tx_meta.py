#!/usr/bin/env python3
import sys, re, gzip
gtf, out = sys.argv[1], sys.argv[2]
g = lambda a, k: (re.search(k + r' "([^"]+)"', a) or [None, None])[1]
ex = {}   # transcript -> exon count
rows = {}
with gzip.open(gtf, "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); a = f[8]
        if f[2] == "transcript":
            t = g(a, "transcript_id")
            rows[t] = [t, t.split(".")[0], g(a, "gene_id").split(".")[0], g(a, "gene_name") or "",
                       g(a, "gene_type") or "", g(a, "transcript_type") or "",
                       "1" if 'tag "MANE_Select"' in a else "0",
                       "1" if 'tag "Ensembl_canonical"' in a else "0",
                       g(a, "transcript_support_level") or "NA", f[0], f[3], f[4], f[6]]
        elif f[2] == "exon":
            ex[g(a, "transcript_id")] = ex.get(g(a, "transcript_id"), 0) + 1
with open(out, "w") as o:
    o.write("\t".join(["transcript_id_v","transcript_id","gene_id","gene_name","gene_type","transcript_type",
                       "mane_select","ensembl_canonical","tsl","hs_chr","hs_start","hs_end","hs_strand","n_exons"]) + "\n")
    for t, r in rows.items(): o.write("\t".join(r + [str(ex.get(t, 0))]) + "\n")
print("transcripts:", len(rows))
