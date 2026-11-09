set -euo pipefail

THREADS=16
MIN_COV=10
MIN_MEAN_DEPTH=20
MIN_GENOME_COV=70

FASTQ_DIR="/data/hugo/fastq/trimmed"
REF="/data/hugo/fastq/parechovirus_complete_genomes.clean.fasta"

RESULTS_DIR="/data/hugo/fastq/trimmed/results_genomes"
TMP_DIR="/data/hugo/fastq/trimmed/tmp_genomes"

mkdir -p "$RESULTS_DIR"
mkdir -p "$TMP_DIR"

echo -e "Sample\tReference\tMappedReads\tGenomeCoverage\tMeanDepth\tStatus" \
> "$RESULTS_DIR/summary.tsv"

cd "$FASTQ_DIR"

for R1 in *_1.trimmed.fastq *_1.trimmed.fastq.gz
do
    [ -f "$R1" ] || continue

    SAMPLE=$(basename "$R1" _1.trimmed.fastq)
    SAMPLE=$(basename "$SAMPLE" _1.trimmed.fastq.gz)

    R2=${R1/_1.trimmed/_2.trimmed}

    [ -f "$R2" ] || continue

    echo ""
    echo "===================================="
    echo "PROCESSING: $SAMPLE"
    echo "===================================="

   
    # 1. SCREEN AGAINST ALL COMPLETE PARECHOVIRUS GENOMES
 

    minimap2 \
        -t "$THREADS" \
        -ax sr \
        "$REF" \
        "$R1" \
        "$R2" 2>/dev/null |
    samtools view -F 4 |
    cut -f3 |
    sort |
    uniq -c |
    sort -nr \
    > "$TMP_DIR/${SAMPLE}.hits.txt"

    BEST=$(head -1 "$TMP_DIR/${SAMPLE}.hits.txt" | awk '{print $2}')

    if [ -z "$BEST" ]; then
        echo "No parechovirus genome detected"
        echo -e "${SAMPLE}\tNA\t0\t0\t0\tFAIL" \
        >> "$RESULTS_DIR/summary.tsv"
        continue
    fi

    echo "Best reference: $BEST"
  
    # 2. EXTRACT BEST COMPLETE GENOME

    seqkit grep \
        -p "$BEST" \
        "$REF" \
        > "$TMP_DIR/${SAMPLE}.ref.fasta"

    samtools faidx "$TMP_DIR/${SAMPLE}.ref.fasta"

    bwa index "$TMP_DIR/${SAMPLE}.ref.fasta" >/dev/null 2>&1


    # 3. REMAP AGAINST BEST COMPLETE GENOME

    bwa mem \
        -t "$THREADS" \
        "$TMP_DIR/${SAMPLE}.ref.fasta" \
        "$R1" \
        "$R2" |
    samtools sort \
        -@ "$THREADS" \
        -o "$TMP_DIR/${SAMPLE}.bam"

    samtools index "$TMP_DIR/${SAMPLE}.bam"

    # 4. COVERAGE STATS

    samtools coverage \
        "$TMP_DIR/${SAMPLE}.bam" \
        > "$TMP_DIR/${SAMPLE}.coverage.txt"

    READS=$(tail -1 "$TMP_DIR/${SAMPLE}.coverage.txt" | awk '{print $4}')
    GENOME_COV=$(tail -1 "$TMP_DIR/${SAMPLE}.coverage.txt" | awk '{print $6}')
    MEAN_DEPTH=$(tail -1 "$TMP_DIR/${SAMPLE}.coverage.txt" | awk '{print $7}')

    samtools depth -aa "$TMP_DIR/${SAMPLE}.bam" \
    > "$TMP_DIR/${SAMPLE}.depth.txt"

    echo "Mapped reads    : $READS"
    echo "Genome coverage : $GENOME_COV"
    echo "Mean depth      : $MEAN_DEPTH"

    # 5. FILTER

    PASS=$(awk -v d="$MEAN_DEPTH" -v c="$GENOME_COV" \
    -v min_d="$MIN_MEAN_DEPTH" -v min_c="$MIN_GENOME_COV" '
    BEGIN{
        if(d>=min_d && c>=min_c)
            print 1;
        else
            print 0;
    }')

    if [ "$PASS" -eq 0 ]; then
        echo "FAILED"
        echo -e "${SAMPLE}\t${BEST}\t${READS}\t${GENOME_COV}\t${MEAN_DEPTH}\tFAIL" \
        >> "$RESULTS_DIR/summary.tsv"
        continue
    fi

    echo "PASSED"
    # 6. MASK LOW-COVERAGE REGIONS
    # Positions with coverage < MIN_COV become N

    awk -v c="$MIN_COV" \
    '$3<c {print $1"\t"$2-1"\t"$2}' \
    "$TMP_DIR/${SAMPLE}.depth.txt" \
    > "$TMP_DIR/${SAMPLE}.lowcov.bed"

  
    # 7. VARIANT CALLING
  
    bcftools mpileup \
        -Ou \
        -f "$TMP_DIR/${SAMPLE}.ref.fasta" \
        -Q 20 \
        "$TMP_DIR/${SAMPLE}.bam" |
    bcftools call \
        -mv \
        -Oz \
        -o "$TMP_DIR/${SAMPLE}.vcf.gz"

    bcftools index "$TMP_DIR/${SAMPLE}.vcf.gz"


    # 8. CONSENSUS GENOME
    # Low-coverage regions are masked as N

    bcftools consensus \
        -f "$TMP_DIR/${SAMPLE}.ref.fasta" \
        -m "$TMP_DIR/${SAMPLE}.lowcov.bed" \
        "$TMP_DIR/${SAMPLE}.vcf.gz" \
        > "$TMP_DIR/${SAMPLE}.consensus.tmp.fasta"

    # 9. HEADER RENAMING
    
    REF_HEADER=$(grep "^>" "$TMP_DIR/${SAMPLE}.ref.fasta" | sed 's/^>//')

    VIRUS=$(echo "$REF_HEADER" | cut -d'|' -f2)
    GENOTYPE=$(echo "$REF_HEADER" | cut -d'|' -f3)

    {
        echo ">${SAMPLE}|${VIRUS}|${GENOTYPE}|complete_genome"
        grep -v "^>" "$TMP_DIR/${SAMPLE}.consensus.tmp.fasta"
    } > "$RESULTS_DIR/${SAMPLE}.complete_genome.fasta"

    rm "$TMP_DIR/${SAMPLE}.consensus.tmp.fasta"

    echo -e "${SAMPLE}\t${BEST}\t${READS}\t${GENOME_COV}\t${MEAN_DEPTH}\tPASS" \
    >> "$RESULTS_DIR/summary.tsv"

done
