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
# import anndata as ad

import pickle

# ─────────────────────────────────────────────
# 0. CONFIGURATION — edit these paths/params
# ─────────────────────────────────────────────
SAMPLE_NAME     = "t2"
IN_DIR          = f"/home/qlyu/mydata/data/v1_astro/{SAMPLE_NAME}/spin"
WORK_DIR        = f"/home/qlyu/mydata/data/v1_astro/{SAMPLE_NAME}/spout"
DB_DIR          = "/home/qlyu/mydata/data/common"
N_TOPICS        = [10, 20, 30]   # LDA topic range for pycisTopic
N_CPUS          = 8

CELLTYPE_COL    = "Subclass"
PATH_RNA_H5AD   = f"{IN_DIR}/rna_preprocessed_astro_yoo25_v2.h5ad"    
PATH_ATAC_FRAG  = f"{IN_DIR}/merged_fragments_v2.tsv.gz"          
PATH_REGIONS    = f"{IN_DIR}/t2_summits_501_filtered.bed" 
PATH_BLACKLIST  = f"{DB_DIR}/mm10.blacklist.bed"
PATH_OUT_REGION_SETS = f'{IN_DIR}/region_sets'

os.makedirs(WORK_DIR, exist_ok=True)

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import find_diff_features, impute_accessibility


# # ─────────────────────────────────────────────
# # 2. LOAD AND PREPROCESS ATAC  (pycisTopic)
# # ─────────────────────────────────────────────
# print("\n=== Step 2: Build pycisTopic object from ATAC ===")


# cistopic_obj = create_cistopic_object_from_fragments(
#     path_to_fragments=PATH_ATAC_FRAG,
#     path_to_regions=PATH_REGIONS,          
#     path_to_blacklist=PATH_BLACKLIST,        
#     project=SAMPLE_NAME,
#     # split_pattern="-", # THIS will create duplicate cells - don't do this
# )
# print(f"  ATAC: {cistopic_obj.binary_matrix.shape}")

# with open(os.path.join(WORK_DIR, "cistopic_obj.pkl"), "wb") as f:
#     pickle.dump(cistopic_obj, f)

# # ─────────────────────────────────────────────
# # 3. RUN LDA TOPIC MODELLING (pycisTopic)
# # ─────────────────────────────────────────────
# print("\n=== Step 3: LDA topic modelling ===")

# with open(os.path.join(WORK_DIR, "cistopic_obj.pkl"), "rb") as f:
#     cistopic_obj = pickle.load(f)

# models = run_cgs_models(
#     cistopic_obj,
#     n_topics=N_TOPICS,
#     n_cpu=N_CPUS,
#     n_iter=10,
#     random_state=555,
#     alpha=50,
#     alpha_by_topic=True,
#     eta=0.1,
#     eta_by_topic=False,
#     save_path=os.path.join(WORK_DIR, "lda_models"),
# )

# # Select best model (highest coherence)
# from pycisTopic.lda_models import evaluate_models
# model = evaluate_models(
#     models,
#     select_model=None,   # None → auto-select
#     return_model=True,
#     metrics=["Minmo_2011", "loglikelihood"],
#     plot=True,
#     save=os.path.join(WORK_DIR, "model_evaluation.pdf"),
# )
# cistopic_obj.add_LDA_model(model)

# with open(os.path.join(WORK_DIR, "cistopic_obj_model.pkl"), "wb") as f:
#     pickle.dump(cistopic_obj, f)

with open(os.path.join(WORK_DIR, "cistopic_obj_model.pkl"), "rb") as f:
    cistopic_obj = pickle.load(f)


# Binarize topics → get topic-specific open regions
region_bin_topics_top3k = binarize_topics(
    cistopic_obj, method="ntop", ntop=3000, plot=True,
    save=os.path.join(WORK_DIR, "binarized_topics_top3k.pdf"),
)

# # 2. Save each topic as a separate BED file
# TOPICS_SUBFOLDER = 'topics_top_3k'

# os.makedirs(PATH_OUT_REGION_SETS, exist_ok=True)
# os.makedirs(os.path.join(PATH_OUT_REGION_SETS, TOPICS_SUBFOLDER), exist_ok=True)
# for topic, regions in region_bin_topics_top3k.items():
#     # regions is a list of strings like "chr1:100-200"
#     regions = regions.index.values
#     with open(f'{PATH_OUT_REGION_SETS}/{TOPICS_SUBFOLDER}/{topic}.bed', 'w') as f:
#         for r in regions:
#             chrom, rest = r.split(':')
#             start, end = rest.split('-')
#             f.write(f"{chrom}\t{start}\t{end}\n")


# ─────────────────────────────────────────────
# 4. DIFFERENTIAL ACCESSIBLE REGIONS
# ─────────────────────────────────────────────
print("\n=== Step 4: Differential accessible regions per cell type ===")

rna = sc.read(PATH_RNA_H5AD, backed='r') # only to read obs
# Transfer cluster labels from RNA to ATAC
cell_data = rna.obs[[CELLTYPE_COL]].copy()
cell_data.index = cell_data.index  # must match ATAC barcodes
cistopic_obj.add_cell_data(cell_data)

imputed_features_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6,
    chunk_size=20000,
    project="cisTopic_Impute",
)

dar_dict = find_diff_features(
    cistopic_obj,
    imputed_features_obj,
    variable=CELLTYPE_COL,
    var_features=None, # region_bin_topics_top3k,
    contrasts=None,   # None → all pairwise
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=N_CPUS,
)

# 2. Save each topic as a separate BED file
TOPICS_SUBFOLDER = 'topics_dar'

os.makedirs(PATH_OUT_REGION_SETS, exist_ok=True)
os.makedirs(os.path.join(PATH_OUT_REGION_SETS, TOPICS_SUBFOLDER), exist_ok=True)
for topic, regions in dar_dict.items():
    # regions is a list of strings like "chr1:100-200"
    regions = regions.index.values
    with open(f'{PATH_OUT_REGION_SETS}/{TOPICS_SUBFOLDER}/{topic.replace("/","")}.bed', 'w') as f:
        for r in regions:
            chrom, rest = r.split(':')
            start, end = rest.split('-')
            f.write(f"{chrom}\t{start}\t{end}\n")