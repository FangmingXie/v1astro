"""
SCENIC+ Pipeline for Single-Cell Multiome Data (scRNA-seq + scATAC-seq)
========================================================================

Reference: https://scenicplus.readthedocs.io
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ─────────────────────────────────────────────
# 0. CONFIGURATION — edit these paths/params
# ─────────────────────────────────────────────
SAMPLE_NAME     = "t0"
WORK_DIR        = "/home/qlyu/mydata/data/v1_astro/scenicplus_output"
RNA_H5_PATH     = "/home/qlyu/mydata/data/v1_astro/scenicplus_inputdata/rna_preprocessed_astro_yoo25.h5ad"   # 10x Cell Ranger ARC output
ATAC_FRAG_PATH  = "/home/qlyu/mydata/data/v1_astro/scenicplus_inputdata/astro_merged_fragments.tsv.gz"       # 10x fragment file
GENOME          = "mm10"                             # or "hg38"
CISTROME_DB     = "/home/qlyu/mydata/data/v1_astro/scenicplus_databases/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather"   # path to cistromes DB
N_TOPICS        = [10, 20, 30]   # LDA topic range for pycisTopic
N_CPUS          = 8

PATH_TO_GENOME_ANNOTATION = "/home/qlyu/mydata/data/v1_astro/scenicplus_databases/mm10.chrom.sizes"
PATH_TO_MOTIF_ANNOTATIONS = "/home/qlyu/mydata/data/v1_astro/scenicplus_databases/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

os.makedirs(WORK_DIR, exist_ok=True)

# ─────────────────────────────────────────────
# 1. LOAD AND PREPROCESS RNA
# ─────────────────────────────────────────────
# print("=== Step 1: Load & preprocess scRNA-seq ===")

# # Load 10x multiome RNA data
# rna = sc.read_10x_h5(RNA_H5_PATH, gex_only=True)
# rna.var_names_make_unique()

# # Basic QC
# sc.pp.filter_cells(rna, min_genes=200)
# sc.pp.filter_genes(rna, min_cells=3)
# rna.var["mt"] = rna.var_names.str.startswith("mt-")  # mouse uses lowercase "mt-"
# sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], inplace=True)
# rna = rna[rna.obs.pct_counts_mt < 20].copy()

# # Normalise & cluster
# sc.pp.normalize_total(rna, target_sum=1e4)
# sc.pp.log1p(rna)
# sc.pp.highly_variable_genes(rna, n_top_genes=3000)
# rna.raw = rna
# rna = rna[:, rna.var.highly_variable].copy()
# sc.pp.scale(rna, max_value=10)
# sc.tl.pca(rna, n_comps=50)
# sc.pp.neighbors(rna)
# sc.tl.umap(rna)
# sc.tl.leiden(rna, resolution=0.5)

# rna.write_h5ad(os.path.join(WORK_DIR, "rna_preprocessed.h5ad"))
# print(f"  RNA: {rna.shape[0]} cells × {rna.shape[1]} genes")

# ─────────────────────────────────────────────
# 2. LOAD AND PREPROCESS ATAC  (pycisTopic)
# ─────────────────────────────────────────────
print("\n=== Step 2: Build pycisTopic object from ATAC ===")

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import (
    find_diff_features,
    # get_conservation_scores,
)

rna = sc.read(RNA_H5_PATH)

# Cell barcodes shared with RNA
barcodes = list(rna.obs_names)

# Genome-specific chromosome sizes (needed to create pseudobulk bigwigs)
if GENOME == "mm10":
    chromsizes_url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
elif GENOME == "hg38":
    chromsizes_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=ATAC_FRAG_PATH,
    path_to_regions=None,          # None → call peaks with MACS2 internally
    path_to_blacklist=None,        # optional BED file
    # split_pattern="-", # THIS will create duplicate cells - don't do this
    project=SAMPLE_NAME,
)
print(f"  ATAC: {cistopic_obj.binary_matrix.shape}")

# Filter to cells shared with RNA
shared_barcodes = list(
    set(cistopic_obj.cell_names) & set(barcodes)
)
cistopic_obj.subset(cells=shared_barcodes)
print(f"  Shared barcodes: {len(shared_barcodes)}")

# ─────────────────────────────────────────────
# 3. RUN LDA TOPIC MODELLING (pycisTopic)
# ─────────────────────────────────────────────
print("\n=== Step 3: LDA topic modelling ===")

models = run_cgs_models_mallet(
    cistopic_obj,
    n_topics=N_TOPICS,
    n_cpu=N_CPUS,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(WORK_DIR, "mallet_tmp"),
    save_path=os.path.join(WORK_DIR, "lda_models"),
)

# Select best model (highest coherence)
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model=None,   # None → auto-select
    return_model=True,
    metrics=["Minmo_2011", "loglikelihood"],
    plot=True,
    save=os.path.join(WORK_DIR, "model_evaluation.pdf"),
)
cistopic_obj.add_LDA_model(model)

# Binarize topics → get topic-specific open regions
region_bin_topics_top3k = binarize_topics(
    cistopic_obj, method="top3k", plot=True,
    save=os.path.join(WORK_DIR, "binarized_topics_top3k.pdf"),
)
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method="otsu", plot=True,
    save=os.path.join(WORK_DIR, "binarized_topics_otsu.pdf"),
)

import pickle
with open(os.path.join(WORK_DIR, "cistopic_obj.pkl"), "wb") as f:
    pickle.dump(cistopic_obj, f)

# ─────────────────────────────────────────────
# 4. DIFFERENTIAL ACCESSIBLE REGIONS
# ─────────────────────────────────────────────
print("\n=== Step 4: Differential accessible regions per cell type ===")


# Transfer cluster labels from RNA to ATAC
cell_data = rna.obs[["leiden"]].copy()
cell_data.index = cell_data.index  # must match ATAC barcodes
cistopic_obj.add_cell_data(cell_data)

imputed_acc_obj = find_diff_features(
    cistopic_obj,
    variable="leiden",
    var_features=region_bin_topics_otsu,
    contrasts=None,   # None → all pairwise
    adjpvalue_thr=0.05,
    log2fc_thr=0.5,
    n_cpu=N_CPUS,
)

# ─────────────────────────────────────────────
# 5. MOTIF ENRICHMENT  (pycistarget)
# ─────────────────────────────────────────────
print("\n=== Step 5: Motif enrichment (pycistarget) ===")

from pycistarget.motif_enrichment_dem import dem
from pycistarget.motif_enrichment_cistarget import run_cistarget

# mm10 motif annotation file (downloaded in setup_scenicplus_env.sh):
#   motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
# Point CISTROME_DB to the rankings .feather file:
#   mm10_screen_v10_clust.regions_vs_motifs.rankings.feather
dem_result = dem(
    region_sets=region_bin_topics_top3k,
    path_to_motif_annotations=CISTROME_DB,   # auto-fetched if None
    genome_annotation=GENOME,
    n_cpu=N_CPUS,
    path_tmp=os.path.join(WORK_DIR, "dem_tmp"),
    save_path=os.path.join(WORK_DIR, "dem_results"),
)

cistarget_result = run_cistarget(
    region_sets=region_bin_topics_top3k,
    rankings_db=None,       # auto-fetched if None
    motif_annotation=None,
    genome_annotation=GENOME,
    n_cpu=N_CPUS,
    path_tmp=os.path.join(WORK_DIR, "cistarget_tmp"),
    save_path=os.path.join(WORK_DIR, "cistarget_results"),
)

# ─────────────────────────────────────────────
# 6. BUILD SCENIC+ OBJECT
# ─────────────────────────────────────────────
print("\n=== Step 6: Create SCENIC+ object ===")

from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata=rna.raw.to_adata(),
    cisTopic_obj=cistopic_obj,
    menr={                          # merged motif enrichment results
        "DEM": dem_result,
        "cisTarget": cistarget_result,
    },
    region_to_gene_approach="GBM",  # Gradient Boosting Machine
    multi_ome_mode=True,
    key_to_group_by="leiden",
    nr_cells_per_metacells=10,
    cell_knn_k=5,
)

# ─────────────────────────────────────────────
# 7. INFER eREGULONS
# ─────────────────────────────────────────────
print("\n=== Step 7: Infer eRegulons ===")

from scenicplus.eregulon_enrichment import score_eRegulons
from scenicplus.wrappers.run_scenicplus import run_scenicplus

run_scenicplus(
    scplus_obj=scplus_obj,
    variable=["leiden"],
    work_dir=WORK_DIR,
    biomart_host="http://www.ensembl.org",
    tf_file=None,        # path to TF list; None = all TFs
    save_format="h5ad",
    prefix=SAMPLE_NAME,
    njobs=N_CPUS,
    # Options: GBM, RF (random forest)
    method_mdl="GBM",
    order_pr=False,
)

# ─────────────────────────────────────────────
# 8. SCORE + FILTER eREGULONS
# ─────────────────────────────────────────────
print("\n=== Step 8: Score & filter eRegulons ===")

score_eRegulons(
    scplus_obj,
    ranking_db_fname=None,   # auto-used from cistarget step
    eRegulon_signatures_key="eRegulon_signatures_filtered",
    key_added="eRegulon_AUC",
    enrichment_type="region",
    n_cpu=N_CPUS,
)

apply_std_filtering_to_eRegulons(scplus_obj)

# ─────────────────────────────────────────────
# 9. DOWNSTREAM: UMAP + DOTPLOT
# ─────────────────────────────────────────────
print("\n=== Step 9: Visualisation ===")

from scenicplus.dimensionality_reduction import (
    run_eRegulons_UMAP,
    run_eRegulons_tSNE,
)
from scenicplus.plotting.dotplot import heatmap_dotplot

run_eRegulons_UMAP(
    scplus_obj,
    auc_key="eRegulon_AUC",
    reduction_name="eRegulons_UMAP",
    random_state=555,
)

heatmap_dotplot(
    scplus_obj,
    size_matrix=scplus_obj.uns["eRegulon_AUC"]["Region_based"],
    color_matrix=scplus_obj.uns["eRegulon_AUC"]["Gene_based"],
    scale_size_matrix=True,
    scale_color_matrix=True,
    group_variable="leiden",
    figsize=(40, 20),
    save=os.path.join(WORK_DIR, "eRegulon_dotplot.pdf"),
)

# ─────────────────────────────────────────────
# 10. SAVE FINAL OBJECT
# ─────────────────────────────────────────────
import dill
with open(os.path.join(WORK_DIR, f"{SAMPLE_NAME}_scplus_obj.pkl"), "wb") as f:
    dill.dump(scplus_obj, f)

print(f"\n✓ SCENIC+ pipeline complete. Results saved to: {WORK_DIR}")
print("  Key outputs:")
print(f"    {WORK_DIR}/{SAMPLE_NAME}_scplus_obj.pkl   ← final SCENIC+ object")
print(f"    {WORK_DIR}/eRegulon_dotplot.pdf           ← eRegulon activity dotplot")
print(f"    {WORK_DIR}/rna_preprocessed.h5ad          ← pre-processed RNA AnnData")
