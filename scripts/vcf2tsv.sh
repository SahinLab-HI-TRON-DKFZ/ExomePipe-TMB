#!/bin/bash

# Usage:
# bash vcf2tsv.sh <input_vcf> <output_tsv>

set -e

input_vcf="$1"
output_tsv="$2"

# === Validation ===
if [[ ! -f "$input_vcf" ]]; then
    echo "[ERROR] Input VCF not found: $input_vcf"
    exit 1
fi

echo "[INFO] Converting VCF to TSV..."
vcf2tsv "$input_vcf" > "$output_tsv"

echo "[INFO] Conversion complete: $input_vcf → $output_tsv"
