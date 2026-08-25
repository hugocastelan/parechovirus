#!/bin/bash

OUTDIR="/data/hugo/parechovirus_missing_fastq"
mkdir -p "$OUTDIR"

while read SRA
do
    echo "Downloading $SRA..."

    prefetch "$SRA" --output-directory "$OUTDIR/sra"

    fasterq-dump \
        "$OUTDIR/sra/$SRA/$SRA.sra" \
        --split-files \
        --threads 8 \
        --outdir "$OUTDIR"

    echo "$SRA finished"

done < /home/hugocastelan/Documents/projects/human_parechovirus/new_id_Serratus_database/id_misssing_sra_parechovirus_all.txt
