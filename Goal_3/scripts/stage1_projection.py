#!/usr/bin/env python3
"""Human gene -> chimp gene projection table. One mutually exclusive anchor category per human gene."""
import argparse, gzip, re, collections
import pandas as pd
p = argparse.ArgumentParser()
for a in ["gencode-gtf","entrez-meta","gene-info","orthologs","chimp-genes","in-annot","out"]: p.add_argument("--"+a, required=True)
A = p.parse_args()
strip = lambda s: re.sub(r"\.\d+$", "", s)

genes, tx2gene = {}, {}
with gzip.open(A.gencode_gtf, "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); a = f[8]
        gid = strip(re.search(r'gene_id "([^"]+)"', a).group(1))
        if f[2] == "gene":
            genes[gid] = dict(human_gene_id=gid, human_gene_symbol=(re.search(r'gene_name "([^"]+)"', a) or [0,""])[1],
                              human_gene_biotype=(re.search(r'gene_type "([^"]+)"', a) or [0,""])[1])
        elif f[2] == "transcript":
            tx2gene[strip(re.search(r'transcript_id "([^"]+)"', a).group(1))] = gid
H = pd.DataFrame(genes.values())

g2e = collections.defaultdict(set)                      # source 1: GENCODE metadata (transcript level -> gene)
with gzip.open(A.entrez_meta, "rt") as fh:
    for l in fh:
        t, e = l.rstrip("\n").split("\t")[:2]
        g = tx2gene.get(strip(t))
        if g: g2e[g].add(e)
g2e_ncbi, sym2id = collections.defaultdict(set), collections.defaultdict(set)   # source 2: NCBI gene_info Ensembl xref
with gzip.open(A.gene_info, "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); sym2id[f[2]].add(f[1])
        for m in re.findall(r"Ensembl:(ENSG\d+)", f[5]): g2e_ncbi[m].add(f[1])

def link(r):
    s1, s2 = g2e.get(r.human_gene_id, set()), g2e_ncbi.get(r.human_gene_id, set())
    both = s1 | s2
    if len(both) == 1: return next(iter(both)), "gencode+ncbi" if (s1 and s2) else ("gencode" if s1 else "ncbi_xref")
    if len(both) > 1:
        ag = s1 & s2
        if len(ag) == 1: return next(iter(ag)), "agree_of_multiple"
        return "", "ambiguous:" + ";".join(sorted(both))
    s = sym2id.get(r.human_gene_symbol, set())
    if len(s) == 1: return next(iter(s)), "symbol_rescue"
    return "", "none"
H[["human_geneid","geneid_link_source"]] = H.apply(lambda r: pd.Series(link(r)), axis=1)

O = pd.read_csv(A.orthologs, sep="\t", dtype=str); O.columns = [c.lstrip("#") for c in O.columns]
O = O[(O.tax_id == "9606") & (O.Other_tax_id == "9598")]
h2c = O.groupby("GeneID").Other_GeneID.apply(lambda s: sorted(set(s))).to_dict()
c2h = O.groupby("Other_GeneID").GeneID.apply(lambda s: sorted(set(s))).to_dict()
C = pd.read_csv(A.chimp_genes, sep="\t", dtype=str)
ndup = int(C.chimp_geneid.duplicated().sum())
if ndup:
    C["span"] = C.chimp_end.astype(int) - C.chimp_start.astype(int)
    C = C.sort_values("span", ascending=False).drop_duplicates("chimp_geneid")
    print(f"WARNING: {ndup} duplicate GeneID gene records in chimp GTF; kept the longest record each")
C = C.set_index("chimp_geneid")
in_annot = set(l.strip() for l in open(A.in_annot))

def classify(r):
    if not r.human_geneid: return "mapping_failure", [], ""
    cg = h2c.get(r.human_geneid, [])
    if not cg: return "no_ncbi_ortholog", [], ""
    known = [c for c in cg if c in C.index]
    if not known: return "ortholog_not_in_refseq_gtf", [], "GeneID(s) not in chimp GTF: " + ";".join(cg)
    if len(known) > 1: return "one_to_many", known, ""
    others = [h for h in c2h.get(known[0], []) if h != r.human_geneid]
    if others: return "many_to_one_complex", known, "other human GeneIDs on same chimp gene: " + ";".join(others)
    if known[0] not in in_annot: return "ortholog_absent_from_annotation", known, ""
    return "one_to_one", known, ""
H[["orthology_class","chimp_geneids","notes"]] = H.apply(lambda r: pd.Series(classify(r)), axis=1)
H["chimp_gene_id"] = H.chimp_geneids.apply(";".join)
def loci(ids):
    return ";".join(f"{C.loc[c].chimp_chr}:{C.loc[c].chimp_start}-{C.loc[c].chimp_end}({C.loc[c].chimp_strand})"
                    for c in ids if c in C.index)
H["chimp_loci"] = H.chimp_geneids.apply(loci)
H["chimp_gene_symbol"] = H.chimp_geneids.apply(lambda ids: ";".join(C.loc[c].chimp_gene_symbol for c in ids if c in C.index))
conf = {"one_to_one":"high", "one_to_many":"medium", "many_to_one_complex":"medium", "ortholog_absent_from_annotation":"high_anchor_no_annotation"}
H["projection_confidence"] = H.orthology_class.map(conf).fillna("none")
H["orthology_source"] = H.orthology_class.apply(lambda c: "ncbi_ortholog" if c in conf else "none")
H["analysis_set"] = H.orthology_class.isin(conf.keys())
H.drop(columns="chimp_geneids").to_csv(A.out, sep="\t", index=False)
print(H.orthology_class.value_counts(), "\n\nprotein-coding only:\n", H[H.human_gene_biotype=="protein_coding"].orthology_class.value_counts())
