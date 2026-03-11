#!/bin/bash

# ===============================
# CONFIGURACIÓN
# ===============================

INPUT_DIR="/data/hugo"
OUTPUT_DIR="/data/hugo/fastp_trimmed"
THREADS=8

mkdir -p $OUTPUT_DIR

echo "Starting fastp trimming"
echo "Date: $(date)"
echo "================================"

for R1 in ${INPUT_DIR}/*_1.fastq
do
    BASENAME=$(basename $R1 _1.fastq)
    R2=${INPUT_DIR}/${BASENAME}_2.fastq

    if [ ! -f "$R2" ]; then
        echo "WARNING: Missing pair for $BASENAME"
        continue
    fi

    echo "--------------------------------"
    echo "Processing $BASENAME"
    echo "Date: $(date)"

    fastp \
        -i $R1 \
        -I $R2 \
        -o ${OUTPUT_DIR}/${BASENAME}_1.trimmed.fastq \
        -O ${OUTPUT_DIR}/${BASENAME}_2.trimmed.fastq \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --n_base_limit 5 \
        --low_complexity_filter \
        --thread $THREADS \
        --html ${OUTPUT_DIR}/${BASENAME}.html \
        --json ${OUTPUT_DIR}/${BASENAME}.json

done

echo "================================"
echo "fastp finished at: $(date)"
