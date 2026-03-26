#!/bin/bash

SCRIPT_DIR="/home/qlyu/mydata/code/trial/create_cisTarget_databases"
DB_DIR="/home/qlyu/mydata/data/common"
WORK_DIR="/home/qlyu/mydata/data/v1_astro/t2/spin"

GENOME_FASTA="${DB_DIR}/mm10.fa"
CHROMSIZES="${DB_DIR}/mm10.chrom.sizes"
MOTIF_DIR="${DB_DIR}/aertslab_motif_collection/v10nr_clust_public/singletons"

REGION_BED="${WORK_DIR}/t2_summits_501_filtered.bed"

OUT_MOTIF_LIST="${DB_DIR}/motifs.txt"
OUT_FASTA="${WORK_DIR}/t2_summits_501_filtered_1kb_bg_padding.fa"
OUT_CISTARGET_DB="${WORK_DIR}/t2"

ls ${MOTIF_DIR} > ${OUT_MOTIF_LIST}
echo "saved ${OUT_MOTIF_LIST}"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${OUT_FASTA} \
        1000 \
        yes

${SCRIPT_DIR}/create_cistarget_motif_databases.py \
    -f ${OUT_FASTA} \
    -M ${MOTIF_DIR} \
    -m ${OUT_MOTIF_LIST} \
    -o ${OUT_CISTARGET_DB} \
    --bgpadding 1000 \
    -t 20