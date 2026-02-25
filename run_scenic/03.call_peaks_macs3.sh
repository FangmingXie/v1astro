#!/usr/bin/env bash
# =============================================================================
# MACS3 Peak Calling Pipeline — ATAC-seq | 10x fragments.tsv.gz | mm10
# Usage: bash macs3_atac_pipeline.sh <fragments.tsv.gz> <sample_name> <outdir>
# =============================================================================

set -euo pipefail

BLACKLIST="/home/qlyu/mydata/data/v1_astro/scenicplus_databases/mm10.blacklist.bed"

# --- Arguments ---------------------------------------------------------------
FRAGMENTS="/home/qlyu/mydata/data/v1_astro/scenicplus_inputdata/astro_merged_fragments.tsv.gz" # "${1:-fragments.tsv.gz}"
SAMPLE="t0" # "${2:-my_sample}"
OUTDIR="/home/qlyu/mydata/data/v1_astro/scenicplus_inputdata/macs3_peaks" #"${3:-./macs3_peaks}"

# --- Setup -------------------------------------------------------------------
mkdir -p "${OUTDIR}"

echo "[$(date)] Starting MACS3 ATAC-seq pipeline"
echo "  Fragments : ${FRAGMENTS}"
echo "  Sample    : ${SAMPLE}"
echo "  Output dir: ${OUTDIR}"

# # =============================================================================
# # STEP 1 — Convert fragments.tsv.gz to BED
# # =============================================================================
# echo ""
# echo "[$(date)] Step 1: Preparing BED file from fragments..."

BED="${OUTDIR}/${SAMPLE}_fragments.bed"

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
#   --outdir "${OUTDIR}" \
#   --nomodel \
#   --shift  -75 \
#   --extsize 150 \
#   -q 0.01 \
#   --call-summits \
#   2> "${OUTDIR}/macs3.log"

# echo "  Peaks called: $(wc -l < "${OUTDIR}/${SAMPLE}_peaks.narrowPeak")"

# =============================================================================
# STEP 2.5 - Extend Summit
# =============================================================================

awk 'BEGIN{OFS="\t"} /^chr/ {$2 = $2 - 250; $3 = $3 + 250; print $1, $2, $3, $4 }' "${OUTDIR}/${SAMPLE}_summits.bed" > "${OUTDIR}/${SAMPLE}_summits_501.bed"


# =============================================================================
# STEP 3 — Download mm10 blacklist and filter peaks
# =============================================================================
echo ""
echo "[$(date)] Step 3: Filtering blacklist regions..."

BEFORE="${OUTDIR}/${SAMPLE}_summits_501.bed"
FILTERED="${OUTDIR}/${SAMPLE}_summits_501_filtered.bed"

bedtools intersect \
  -v \
  -a "${BEFORE}" \
  -b "${BLACKLIST}" \
  > "${FILTERED}"

echo "  Peaks before filtering : $(wc -l < "${BEFORE}")"
echo "  Peaks after filtering  : $(wc -l < "${FILTERED}")"

# # =============================================================================
# # STEP 4 — Convert bedgraph to bigwig (optional, requires bedGraphToBigWig)
# # =============================================================================
# echo ""
# echo "[$(date)] Step 4: Converting pileup bedgraph to bigwig..."

# CHROMSIZES="${OUTDIR}/mm10.chrom.sizes"
# PILEUP="${OUTDIR}/${SAMPLE}_treat_pileup.bdg"
# BIGWIG="${OUTDIR}/${SAMPLE}.bw"

# if command -v fetchChromSizes &> /dev/null && command -v bedGraphToBigWig &> /dev/null; then
#   fetchChromSizes mm10 > "${CHROMSIZES}"
#   bedtools sort -i "${PILEUP}" \
#     | bedGraphToBigWig - "${CHROMSIZES}" "${BIGWIG}"
#   echo "  BigWig written: ${BIGWIG}"
# else
#   echo "  Skipping bigwig conversion (fetchChromSizes or bedGraphToBigWig not found)"
# fi

# =============================================================================
# STEP 5 — Summary
# =============================================================================
echo ""
echo "[$(date)] Pipeline complete. Output files:"
echo ""
echo "  Main peaks (filtered)   : ${FILTERED}"
echo "  Raw peaks               : ${OUTDIR}/${SAMPLE}_peaks.narrowPeak"
echo "  Summits                 : ${OUTDIR}/${SAMPLE}_summits.bed"
# echo "  Signal track (bedgraph) : ${PILEUP}"
# [ -f "${BIGWIG}" ] && echo "  Signal track (bigwig)   : ${BIGWIG}"
echo "  MACS3 log               : ${OUTDIR}/macs3.log"
echo ""
echo "narrowPeak columns: chr | start | end | name | score | strand |"
echo "                    fold_change | -log10(pval) | -log10(qval) | summit_offset"
