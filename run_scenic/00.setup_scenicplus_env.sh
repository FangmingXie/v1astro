#!/usr/bin/env bash
# =============================================================================
# SCENIC+ Environment Setup — mm10 (Mouse)
# =============================================================================
# SCENIC+ is NOT on PyPI. It must be cloned from GitHub and built locally.
# Run this script section-by-section; don't pipe it blindly to bash.
#
# Tested on: Linux (Ubuntu 20.04 / Rocky 8), Python 3.11.0, conda ≥ 23.x
# Reference: https://scenicplus.readthedocs.io/en/latest/install.html
#            https://github.com/aertslab/scenicplus
# =============================================================================

set -euo pipefail

# Verify
python -c "import scenicplus; print('scenicplus:', scenicplus.__version__)"
python -c "import pycisTopic; print('pycisTopic OK')"
python -c "import pycistarget; print('pycistarget OK')"
python -c "import pybedtools; print('pybedtools OK')"

# # ── 5. Install MALLET (required for LDA topic modelling in pycisTopic) ────────
# conda install -c conda-forge mallet -y

# # ── 6. Download mm10 databases ────────────────────────────────────────────────
# # These files are large (several GB). Download to a dedicated directory.
# mkdir -p scenicplus_databases && cd scenicplus_databases

# # cisTarget motif rankings DB — mm10
# wget -c https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather

# # cisTarget motif scores DB — mm10 (optional but recommended)
# wget -c https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather

# # Motif-to-TF annotation table — mouse MGI gene symbols
# wget -c https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

# # mm10 chromosome sizes (used by pycisTopic for bigwig generation)
# wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# cd ..

# echo ""
# echo "=== Database files downloaded to: scenicplus_databases/ ==="
# ls -lh scenicplus_databases/

# # ── 7. Update paths in the pipeline script ────────────────────────────────────
# # Edit scenicplus_multiome_pipeline.py and set:
# #
# #   CISTROME_DB  = "scenicplus_databases/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather"
# #   GENOME       = "mm10"
# #
# # Also update the dem() and run_cistarget() calls to pass the annotation file:
# #   path_to_motif_annotations = "scenicplus_databases/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# echo ""
# echo "=== Setup complete. Activate environment and run pipeline: ==="
# echo "    conda activate scenicplus"
# echo "    python scenicplus_multiome_pipeline.py"

# # ── 8. Expected input files (10x Cell Ranger ARC output) ──────────────────────
# #   filtered_feature_bc_matrix.h5    ← joint RNA+ATAC barcode-feature matrix
# #   atac_fragments.tsv.gz            ← ATAC fragment file
# #   atac_fragments.tsv.gz.tbi        ← tabix index (must exist alongside .gz)