#!/usr/bin/env python3
import sys, re, pandas as pd
proj, meta, fai, out, nchunks = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
FLANK = 10000
P = pd.read_csv(proj, sep="\t", dtype=str).fillna("")
P = P[P.analysis_set == "True"]
M = pd.read_csv(meta, sep="\t", dtype=str)
clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
tx_by_gene = M.groupby("gene_id").transcript_id.apply(list).to_dict()
rows = []
for _, r in P.iterrows():
    txs = tx_by_gene.get(r.human_gene_id, [])
    if not txs: continue
    loci = []
    for cg, loc in zip(r.chimp_gene_id.split(";"), r.chimp_loci.split(";")):
        m = re.match(r"(.+):(\d+)-(\d+)\(([+-.])\)", loc)
        ch, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        loci.append(f"{cg}|{ch}|{max(1, s-FLANK)}|{min(clen[ch], e+FLANK)}|{m.group(4)}|{s}|{e}")
    rows.append([r.human_gene_id, r.orthology_class, ",".join(txs), ";".join(loci)])
J = pd.DataFrame(rows, columns=["human_gene_id","orthology_class","human_transcripts","loci"])
J["chunk"] = [i % nchunks for i in range(len(J))]
J.to_csv(out, sep="\t", index=False)
print("genes:", len(J), " transcripts:", sum(len(t.split(",")) for t in J.human_transcripts), " chunks:", nchunks)
