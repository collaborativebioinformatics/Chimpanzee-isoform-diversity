#!/usr/bin/env python3
import sys, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
corr_f, proj_f, outd = sys.argv[1:4]
D = pd.read_csv(corr_f, sep="\t", dtype=str); P = pd.read_csv(proj_f, sep="\t", dtype=str)
for c in ["qcov","identity","mapq","shared_junctions","n_exons","n_exons_proj","tss_difference_bp","tes_difference_bp"]:
    if c in D: D[c] = pd.to_numeric(D[c], errors="coerce")
order = ["exact_intron_chain","same_chain_variable_ends","compatible_contained",
         "not_annotated_known_junctions","not_annotated_novel_junctions","divergent_structure",
         "no_annotated_chimp_isoform","locus_aligned_ambiguous","unalignable_in_locus","no_gene_anchor"]

# 1. gene-anchor funnel
pc = P[P.human_gene_biotype == "protein_coding"]
funnel = pd.DataFrame({"all_genes":[len(P), len(pc)],
    "with_GeneID":[(P.human_geneid.notna()&(P.human_geneid!="")).sum(), (pc.human_geneid.notna()&(pc.human_geneid!="")).sum()],
    "ncbi_1to1":[(P.orthology_class=="one_to_one").sum(), (pc.orthology_class=="one_to_one").sum()],
    "1many_or_complex":[P.orthology_class.isin(["one_to_many","many_to_one_complex"]).sum(), pc.orthology_class.isin(["one_to_many","many_to_one_complex"]).sum()],
    "unanchored":[(~P.orthology_class.isin(["one_to_one","one_to_many","many_to_one_complex","ortholog_absent_from_annotation"])).sum(),
                  (~pc.orthology_class.isin(["one_to_one","one_to_many","many_to_one_complex","ortholog_absent_from_annotation"])).sum()]},
    index=["all_biotypes","protein_coding"]).T
funnel.to_csv(f"{outd}/gene_anchor_funnel.tsv", sep="\t"); funnel.plot(kind="bar", figsize=(8,4)); plt.ylabel("human genes"); plt.tight_layout(); plt.savefig(f"{outd}/gene_anchor_funnel.png", dpi=150); plt.close()

# 2. transcript waterfall
w = D.splice_match_category.value_counts().reindex(order, fill_value=0)
w.to_csv(f"{outd}/transcript_waterfall.tsv", sep="\t"); w.plot(kind="barh", figsize=(8,4)); plt.xlabel("human transcripts"); plt.tight_layout(); plt.savefig(f"{outd}/transcript_waterfall.png", dpi=150); plt.close()

# 3. stacked bars by biotype / MANE / orthology class
def stacked(by, name):
    t = D.groupby([by,"splice_match_category"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
    t = t[t.sum(axis=1) > 200]
    pct = t.div(t.sum(axis=1), axis=0).mul(100)
    t.to_csv(f"{outd}/match_by_{name}_counts.tsv", sep="\t"); pct.round(1).to_csv(f"{outd}/match_by_{name}_pct.tsv", sep="\t")
    pct.plot(kind="barh", stacked=True, figsize=(10,5), colormap="tab20"); plt.xlabel("% of human transcripts"); plt.legend(bbox_to_anchor=(1,1), fontsize=7); plt.tight_layout(); plt.savefig(f"{outd}/match_by_{name}.png", dpi=150); plt.close()
D["biotype_coarse"] = D.human_transcript_biotype.apply(lambda t: t if t in ("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay","processed_transcript") else ("pseudogene" if "pseudogene" in str(t) else "other"))
stacked("biotype_coarse","biotype"); stacked("human_canonical_or_mane_status","mane"); stacked("orthology_class","orthology_class")

# 4. alignment quality distributions
A = D[D.projection_tier.isin(["A","B","C"])]
fig, ax = plt.subplots(1, 3, figsize=(14,4))
sns.histplot(A, x="qcov", bins=50, ax=ax[0]); sns.histplot(A, x="identity", bins=50, ax=ax[1]); sns.histplot(A, x="mapq", bins=30, ax=ax[2])
plt.tight_layout(); plt.savefig(f"{outd}/alignment_quality.png", dpi=150); plt.close()

# 5. junction concordance
if "shared_junctions" in D:
    J = A.dropna(subset=["shared_junctions","n_exons_proj"]).copy(); J["proj_junctions"] = J.n_exons_proj - 1
    plt.figure(figsize=(5,5)); plt.hexbin(J.proj_junctions, J.shared_junctions, gridsize=40, bins="log"); plt.plot([0,60],[0,60],"r--",lw=.8)
    plt.xlabel("projected human junctions"); plt.ylabel("junctions shared with chimp hit"); plt.tight_layout(); plt.savefig(f"{outd}/junction_concordance.png", dpi=150); plt.close()

# 6. chromosome QC
ch = A.groupby("chimp_chr").agg(n=("human_transcript_id","size"), ties=("ambiguity_flag", lambda s: s.str.contains("tie").sum()))
ch.to_csv(f"{outd}/chromosome_qc.tsv", sep="\t")

# headline numbers (denominators encoded in the names; never quote without them)
proj = D[~D.splice_match_category.isin(["no_gene_anchor"])]
mane = D[D.human_canonical_or_mane_status == "MANE_Select"]
match3 = ["exact_intron_chain","same_chain_variable_ends","compatible_contained"]
pcg = P[P.human_gene_biotype == "protein_coding"]
head = pd.Series({
  f"MANE_identical_intron_chain_pct_of_{len(mane)}": round(100*mane.splice_match_category.isin(match3[:2]).mean(),1),
  f"MANE_matched_or_compatible_pct_of_{len(mane)}": round(100*mane.splice_match_category.isin(match3).mean(),1),
  f"projected_matched_or_compatible_pct_of_{len(proj)}": round(100*proj.splice_match_category.isin(match3).mean(),1),
  f"all_tx_matched_or_compatible_pct_of_{len(D)}": round(100*D.splice_match_category.isin(match3).mean(),1),
  f"protein_coding_genes_anchored_pct_of_{len(pcg)}": round(100*pcg.analysis_set.eq("True").mean(),1),
  f"all_genes_anchored_pct_of_{len(P)}": round(100*P.analysis_set.eq("True").mean(),1),
  "tx_unresolved_novel_junctions_n": int((D.splice_match_category=="not_annotated_novel_junctions").sum())})
head.to_csv(f"{outd}/headline_numbers.tsv", sep="\t", header=False); print(head)
