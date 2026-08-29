#!/usr/bin/env python3
"""Merge isomatch classification + projection metrics + Stage 1 + tx metadata into the correspondence table."""
import sys, pandas as pd
cls_f, met_f, meta_f, proj_f, chimp_genes_f, out = sys.argv[1:7]

C = pd.read_csv(cls_f, sep="\t", dtype=str)
C["human_transcript_id"] = C.isoform_id.str.replace(r"^HUMANPROJ_", "", regex=True).str.replace(r"__to__.*$", "", regex=True)
C = C.drop_duplicates("human_transcript_id", keep="first")
C = C.rename(columns={"structural_category":"isom_category","subcategory":"isom_subcategory",
    "ref_tx_id":"chimp_transcript_id","ref_gene_id":"chimp_hit_gene_symbol",
    "diff_to_tss":"tss_difference_bp","diff_to_tes":"tes_difference_bp",
    "matched_junctions":"shared_junctions","matched_exons":"shared_exons"})
C = C[["human_transcript_id","isom_category","isom_subcategory","chimp_transcript_id","chimp_hit_gene_symbol",
       "tss_difference_bp","tes_difference_bp","shared_junctions","shared_exons","all_canonical"]]

M = pd.read_csv(met_f, sep="\t", dtype=str)
T = pd.read_csv(meta_f, sep="\t", dtype=str)
P = pd.read_csv(proj_f, sep="\t", dtype=str)[["human_gene_id","human_gene_symbol","chimp_gene_id","chimp_gene_symbol","orthology_class","orthology_source","projection_confidence"]]
D = (T.merge(P, left_on="gene_id", right_on="human_gene_id", how="left")
      .merge(M.drop(columns=["human_gene_id","orthology_class","chimp_gene_id"]), left_on="transcript_id", right_on="human_transcript_id", how="left")
      .merge(C, on="human_transcript_id", how="left"))

def cls(r):
    if pd.isna(r.orthology_class) or r.orthology_class in ("no_ncbi_ortholog","mapping_failure","ortholog_not_in_refseq_gtf"): return "no_gene_anchor"
    if pd.isna(r.projection_tier) or r.projection_tier == "D": return "unalignable_in_locus"
    if r.projection_tier == "C": return "locus_aligned_ambiguous"
    c = str(r.isom_category)
    if c == "full-splice_match":
        return "exact_intron_chain" if str(r.isom_subcategory) == "reference_match" else "same_chain_variable_ends"
    if c == "incomplete-splice_match": return "compatible_contained"
    if c in ("novel_in_catalog","novel_not_in_catalog"): return "divergent_structure"
    return "no_annotated_chimp_isoform"
D["splice_match_category"] = D.apply(cls, axis=1)

D["ambiguity_flag"] = ""
D.loc[D.isom_category.isin(["antisense","intergenic"]), "ambiguity_flag"] = "off_gene"
exp = D.chimp_gene_symbol.fillna("").str.split(";").str[0]
hit_ok = D.chimp_hit_gene_symbol.isna() | (D.chimp_hit_gene_symbol == exp)
D.loc[~hit_ok, "ambiguity_flag"] = (D.loc[~hit_ok, "ambiguity_flag"] + ";hit_outside_anchor_gene").str.lstrip(";")
D["failure_reason"] = D.splice_match_category.where(D.splice_match_category.isin(
    ["no_gene_anchor","unalignable_in_locus","locus_aligned_ambiguous","no_annotated_chimp_isoform"]), "")
D.loc[D.alignment_status == "invalid_projection_zero_len_exon", "failure_reason"] = "invalid_projection"
D["human_canonical_or_mane_status"] = D.mane_select.map({"1":"MANE_Select"}).fillna(D.ensembl_canonical.map({"1":"Ensembl_canonical"})).fillna("other")

cols = ["human_gene_id","human_gene_symbol","transcript_id_v","transcript_type","human_canonical_or_mane_status",
        "chimp_gene_id","chimp_gene_symbol","chimp_transcript_id","chimp_hit_gene_symbol","orthology_source","orthology_class",
        "projection_tier","chimp_chr","chimp_start","chimp_end","chimp_strand","mapq","qcov","identity","n_exons","n_exons_proj",
        "splice_match_category","isom_category","isom_subcategory","shared_junctions","shared_exons",
        "tss_difference_bp","tes_difference_bp","ambiguity_flag","failure_reason"]
D[[c for c in cols if c in D.columns]].rename(columns={"transcript_id_v":"human_transcript_id","transcript_type":"human_transcript_biotype"}).to_csv(out, sep="\t", index=False)
print(D.splice_match_category.value_counts())
print("\nMANE Select only:")
print(D[D.mane_select=="1"].splice_match_category.value_counts())
print("\nhit outside anchor gene:", (~hit_ok).sum())
