#!/usr/bin/env bash
# =============================================================================
# MACS3 Peak Calling Pipeline — ATAC-seq | 10x fragments.tsv.gz | mm10
# Usage: bash macs3_atac_pipeline.sh <fragments.tsv.gz> <sample_name> <outdir>
# =============================================================================

set -euo pipefail

# --- Arguments ---------------------------------------------------------------

SAMPLE="t2" 
DBDIR="/home/qlyu/mydata/data/common"
WORKDIR="/home/qlyu/mydata/data/v1_astro/${SAMPLE}/spin"

BLACKLIST="${DBDIR}/mm10.blacklist.bed"
FRAGMENTS="${WORKDIR}/merged_fragments_v2.tsv.gz"

# --- Setup -------------------------------------------------------------------
mkdir -p "${WORKDIR}"

echo "[$(date)] Starting MACS3 ATAC-seq pipeline"
echo "  Fragments : ${FRAGMENTS}"
echo "  Sample    : ${SAMPLE}"
echo "  Output dir: ${WORKDIR}"

BED="${WORKDIR}/${SAMPLE}_fragments.bed"
SUMMITS="${WORKDIR}/${SAMPLE}_summits.bed"
EXTENDED="${WORKDIR}/${SAMPLE}_summits_501.bed"
FILTERED="${WORKDIR}/${SAMPLE}_summits_501_filtered.bed"


# # # =============================================================================
# # # STEP 1 — Convert fragments.tsv.gz to BED
# # # =============================================================================
# echo ""
# echo "[$(date)] Step 1: Preparing BED file from fragments..."

# zcat "${FRAGMENTS}" \
#   | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2, $3}' \
#   > "${BED}"

# echo "  Written: ${BED}"
# echo "  Total fragments: $(wc -l < "${BED}")"

# # =============================================================================
# # STEP 2 — Call peaks with MACS3
# # =============================================================================
# echo ""
# echo "[$(date)] Step 2: Running MACS3..."

# macs3 callpeak \
#   -f BED \
#   -t "${BED}" \
#   -g mm \
#   -n "${SAMPLE}" \
#   --outdir "${WORKDIR}" \
#   --nomodel \
#   --shift  -75 \
#   --extsize 150 \
#   -q 0.01 \
#   --call-summits \
#   2> "${WORKDIR}/macs3.log"

# echo "  Peaks called: $(wc -l < "${WORKDIR}/${SAMPLE}_peaks.narrowPeak")"

# =============================================================================
# STEP 3 - Extend Summit
# =============================================================================
echo ""
echo "[$(date)] Step 3: Extend summit and filtering blacklist regions..."


awk 'BEGIN{OFS="\t"} /^chr/ {$2 = $2 - 250; $3 = $3 + 250; print $1, $2, $3, $4 }' \
      "${SUMMITS}" \
    > "${EXTENDED}"

bedtools intersect -v -a "${EXTENDED}" -b "${BLACKLIST}" > "${FILTERED}"

# =============================================================================
# STEP 5 — Summary
# =============================================================================
echo ""
echo "[$(date)] Pipeline complete. Output files:"
echo ""
echo "  Main peaks (filtered)   : ${FILTERED}"
echo "  Raw peaks               : ${WORKDIR}/${SAMPLE}_peaks.narrowPeak"
echo "  Summits                 : ${WORKDIR}/${SAMPLE}_summits.bed"
echo "  MACS3 log               : ${WORKDIR}/macs3.log"

