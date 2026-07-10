#!/usr/bin/env bash
# Build a miniature dev dataset for rapid pipeline iteration.
#
# Selects a small number of loci from the full reference and extracts
# read pairs (from each sample) that map to those loci.  The result is
# a self-contained directory that runs each pipeline stage in minutes
# instead of hours.
#
# Output layout (matches working.R expectations):
#   MINI/fastq/      — per-sample FASTQ pairs
#   MINI/reference/  — mini reference FASTA + BWA index
#
# Usage:
#   bash tools/make_mini_data.sh [FULL_DATA_DIR] [MINI_DATA_DIR]
#
# Defaults:
#   FULL_DATA_DIR = dev-data
#   MINI_DATA_DIR = dev-data-mini
#
# Requirements: bwa, samtools, seqtk (all available in the SORTER2 conda env)
# Run from the package root, e.g.:
#   conda run -n SORTER2 bash tools/make_mini_data.sh

set -euo pipefail

FULL="${1:-dev-data}"
MINI="${2:-dev-data-mini}"

REF="${FULL}/fermaily_perfamily2hybpiper_sorter2.fasta"

# --- Parameters -----------------------------------------------------------
# Locus numbers to include (from Stage 1B diploidclusters output; known good)
LOCI=(1 2 3 4 5 6 7 8 9 10)

# Samples: <name>  (expects <name>_R1.fastq and <name>_R2.fastq in FULL)
DIPLOID_SAMPLES=(Iimura12_cthysanostomum)
HYBRID_SAMPLES=(Iimura18_cthysanostomum Iimura46_cgrande)
ALL_SAMPLES=("${DIPLOID_SAMPLES[@]}" "${HYBRID_SAMPLES[@]}")

# Hard cap on read pairs per sample (after extracting mapped reads).
# 0 = no cap.
MAX_READ_PAIRS=30000
# --------------------------------------------------------------------------

MINI_REF_DIR="${MINI}/reference"
MINI_REF="${MINI_REF_DIR}/fermaily_perfamily2hybpiper_sorter2.fasta"
MINI_FASTQ_DIR="${MINI}/fastq"

echo "=== make_mini_data.sh ==="
echo "Full data : ${FULL}"
echo "Mini data : ${MINI}"
echo "Loci      : ${LOCI[*]}"
echo "Samples   : ${ALL_SAMPLES[*]}"
echo ""

mkdir -p "${MINI_REF_DIR}" "${MINI_FASTQ_DIR}"

# --- 1. Build mini reference -----------------------------------------------
echo "[1/3] Building mini reference (${#LOCI[@]} loci)..."

python3 - "${MINI_REF}" "${REF}" ${LOCI[*]} <<'PYEOF'
import sys, re

out_path, ref_path, *loci_args = sys.argv[1:]
loci = set(int(x) for x in loci_args)
pattern = re.compile(r'^>L(\d+)_')
keep = False
with open(ref_path) as fh, open(out_path, 'w') as out:
    for line in fh:
        if line.startswith('>'):
            m = pattern.match(line)
            keep = m is not None and int(m.group(1)) in loci
        if keep:
            out.write(line)
PYEOF

seqcount=$(grep -c "^>" "${MINI_REF}")
echo "   ${seqcount} sequences written to ${MINI_REF}"

# --- 2. Index mini reference -----------------------------------------------
echo "[2/3] Indexing mini reference..."
bwa index "${MINI_REF}" 2>/dev/null

# --- 3. Extract reads per sample -------------------------------------------
echo "[3/3] Extracting reads per sample..."

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    R1_IN="${FULL}/${SAMPLE}_R1.fastq"
    R2_IN="${FULL}/${SAMPLE}_R2.fastq"
    R1_OUT="${MINI_FASTQ_DIR}/${SAMPLE}_R1.fastq"
    R2_OUT="${MINI_FASTQ_DIR}/${SAMPLE}_R2.fastq"

    if [[ ! -f "${R1_IN}" ]]; then
        echo "   SKIP ${SAMPLE}: ${R1_IN} not found"
        continue
    fi

    echo "   ${SAMPLE}: mapping reads to mini reference..."

    TMP_BAM=$(mktemp /tmp/mini_XXXXXX.bam)
    TMP_IDS=$(mktemp /tmp/mini_XXXXXX.ids)

    # Map paired-end reads; extract names of any read where at least one
    # mate mapped (-F 4 keeps records where this read is mapped).
    bwa mem -t "$(nproc)" "${MINI_REF}" "${R1_IN}" "${R2_IN}" 2>/dev/null \
        | samtools view -F 4 \
        | awk '{print $1}' \
        | sort -u \
        > "${TMP_IDS}"

    IDS_COUNT=$(wc -l < "${TMP_IDS}")
    echo "   ${SAMPLE}: ${IDS_COUNT} read names mapped"

    if [[ ${MAX_READ_PAIRS} -gt 0 && ${IDS_COUNT} -gt ${MAX_READ_PAIRS} ]]; then
        echo "   ${SAMPLE}: capping at ${MAX_READ_PAIRS} pairs"
        shuf -n "${MAX_READ_PAIRS}" "${TMP_IDS}" > "${TMP_IDS}.sub"
        mv "${TMP_IDS}.sub" "${TMP_IDS}"
    fi

    # Extract matching read pairs from original FASTQs using seqtk
    seqtk subseq "${R1_IN}" "${TMP_IDS}" > "${R1_OUT}"
    seqtk subseq "${R2_IN}" "${TMP_IDS}" > "${R2_OUT}"

    PAIR_COUNT=$(grep -c "^@" "${R1_OUT}")
    echo "   ${SAMPLE}: ${PAIR_COUNT} pairs written"

    rm -f "${TMP_BAM}" "${TMP_IDS}"
done

# Liu18_cthysanostomum is a synthetic hybrid whose name starts with 'L'.
# It is a copy of Iimura18 and exists solely to exercise the Stage3 rename
# fix (which only triggers when a hybrid sample name begins with 'L').
cp "${MINI_FASTQ_DIR}/Iimura18_cthysanostomum_R1.fastq" \
   "${MINI_FASTQ_DIR}/Liu18_cthysanostomum_R1.fastq"
cp "${MINI_FASTQ_DIR}/Iimura18_cthysanostomum_R2.fastq" \
   "${MINI_FASTQ_DIR}/Liu18_cthysanostomum_R2.fastq"
echo "   Liu18_cthysanostomum: copied from Iimura18 (synthetic L-prefixed hybrid)"

echo ""
echo "Done. Mini dataset: ${MINI}"
echo ""
echo "Suggested test run:"
echo "  Rscript working.R"
