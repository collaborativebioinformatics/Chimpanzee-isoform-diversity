#!/usr/bin/env python3
"""Locus-restricted splice projection of human transcripts. Run per chunk (Slurm array)."""
import sys, os, re, subprocess, tempfile, pysam, pandas as pd
jobs_f, chunk, hs_fa, pt_fa, outdir = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
os.makedirs(outdir, exist_ok=True)
J = pd.read_csv(jobs_f, sep="\t", dtype=str); J = J[J.chunk == str(chunk)]
HS, PT = pysam.FastaFile(hs_fa), pysam.FastaFile(pt_fa)
hs_names = {n.split(".")[0]: n for n in HS.references}
gtf = open(f"{outdir}/proj_chunk{chunk}.gtf", "w")
met = open(f"{outdir}/metrics_chunk{chunk}.tsv", "w")
met.write("\t".join(["human_transcript_id","human_gene_id","chimp_gene_id","orthology_class","n_candidate_loci",
                     "chimp_chr","chimp_start","chimp_end","chimp_strand","n_exons_proj","mapq","AS","qlen","qcov",
                     "identity","n_indels","clip5","clip3","tie","projection_tier","alignment_status"]) + "\n")

def parse(sam_txt, offset, chrom):
    """return list of alignments: dict with exons (genome coords, 1-based), metrics"""
    out = []
    for l in sam_txt.splitlines():
        if not l or l[0] == "@": continue
        f = l.split("\t"); flag = int(f[1])
        if flag & 0x4 or flag & 0x100 or flag & 0x800: continue
        pos = int(f[3]); cig = re.findall(r"(\d+)([MIDNSHX=])", f[5])
        tags = dict((t.split(":")[0], t.split(":")[2]) for t in f[11:])
        exons, ref, cur, q_aln, ins, dele, cols = [], pos, pos, 0, 0, 0, 0
        clip = [0, 0]; first = True
        for n, op in cig:
            n = int(n)
            if op in "M=X": ref += n; q_aln += n; cols += n; first = False
            elif op == "I": q_aln += n; ins += 1; cols += n; first = False
            elif op == "D": ref += n; dele += 1; cols += n; first = False
            elif op == "N": exons.append((cur, ref - 1)); cur = ref + n; ref += n
            elif op in "SH": clip[0 if first else 1] += n
        exons.append((cur, ref - 1))
        qlen = q_aln + clip[0] + clip[1]
        nm = int(tags.get("NM", 0))
        out.append(dict(exons=[(a + offset - 1, b + offset - 1) for a, b in exons],   # offset = locus start (1-based)
                        strand="-" if flag & 0x10 else "+", mapq=int(f[4]), AS=int(tags.get("AS", 0)),
                        qlen=qlen, qcov=q_aln / qlen if qlen else 0, identity=1 - nm / cols if cols else 0,
                        n_indels=ins + dele, clip5=clip[1] if flag & 0x10 else clip[0], clip3=clip[0] if flag & 0x10 else clip[1]))
    return out

with tempfile.TemporaryDirectory() as td:
    for _, r in J.iterrows():
        txs = [t for t in r.human_transcripts.split(",") if t in hs_names]
        if not txs: continue
        qf = f"{td}/q.fa"
        with open(qf, "w") as q:
            for t in txs: q.write(f">{t}\n{HS.fetch(hs_names[t])}\n")
        loci = [x.split("|") for x in r.loci.split(";")]
        best = {t: [] for t in txs}
        for cg, ch, fs, fe, gs, gstart, gend in loci:
            fs, fe = int(fs), int(fe)
            rf = f"{td}/r.fa"
            with open(rf, "w") as o: o.write(f">{ch}:{fs}-{fe}\n{PT.fetch(ch, fs - 1, fe)}\n")
            sam = subprocess.run(["minimap2", "-ax", "splice:hq", "-uf", "--secondary=no", "-t", "1", rf, qf],
                                 capture_output=True, text=True).stdout
            # group by query
            by_q = {}
            for l in sam.splitlines():
                if l and l[0] != "@": by_q.setdefault(l.split("\t")[0], []).append(l)
            for t in txs:
                for a in parse("\n".join(by_q.get(t, [])), fs, ch):
                    a.update(chimp_gene_id=cg, chr=ch); best[t].append(a)
        for t in txs:
            als = sorted(best[t], key=lambda a: -a["AS"])
            if not als:
                met.write("\t".join(map(str, [t, r.human_gene_id, r.loci.split("|")[0], r.orthology_class, len(loci)] + ["NA"] * 15 + ["D", "no_alignment_in_locus"])) + "\n"); continue
            a = als[0]; tie = len(als) > 1 and als[1]["AS"] == a["AS"]
            if tie: tier = "C"
            elif a["qcov"] >= 0.9 and a["identity"] >= 0.95 and a["mapq"] >= 30: tier = "A"
            elif a["qcov"] >= 0.8: tier = "B"
            else: tier = "D"
            status = "projected" if tier in "AB" else ("ambiguous_tie" if tie else "low_coverage")
            pid = f"HUMANPROJ_{t}__to__{a['chimp_gene_id']}"
            attrs = (f'gene_id "{r.human_gene_id}"; transcript_id "{pid}"; human_transcript_id "{t}"; '
                     f'chimp_gene_id "{a["chimp_gene_id"]}"; orthology_class "{r.orthology_class}"; mapq "{a["mapq"]}"; '
                     f'qcov "{a["qcov"]:.3f}"; identity "{a["identity"]:.4f}"; projection_tier "{tier}"; tie "{int(tie)}";')
            if tier != "D":
                ex = a["exons"]
                gtf.write("\t".join([a["chr"], "humanproj", "transcript", str(ex[0][0]), str(ex[-1][1]), ".", a["strand"], ".", attrs]) + "\n")
                for i, (s, e) in enumerate(ex, 1):
                    gtf.write("\t".join([a["chr"], "humanproj", "exon", str(s), str(e), ".", a["strand"], ".", attrs + f' exon_number "{i}";']) + "\n")
            met.write("\t".join(map(str, [t, r.human_gene_id, a["chimp_gene_id"], r.orthology_class, len(loci), a["chr"], a["exons"][0][0], a["exons"][-1][1],
                     a["strand"], len(a["exons"]), a["mapq"], a["AS"], a["qlen"], f'{a["qcov"]:.3f}', f'{a["identity"]:.4f}', a["n_indels"],
                     a["clip5"], a["clip3"], int(tie), tier, status])) + "\n")
gtf.close(); met.close()
