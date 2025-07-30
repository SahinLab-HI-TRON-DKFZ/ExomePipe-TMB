#!/bin/bash

# Usage:
# bash bcftools_norm.sh <input_vcf.gz> <output_vcf.gz> <reference_fasta> [decompressed_output.vcf]

set -e

# Parse arguments
input_vcf="$1"
output_vcf="$2"
reference="$3"
decompressed_vcf="$4"  # optional

if [[ ! -f "$input_vcf" ]]; then
    echo "[ERROR] Input VCF not found: $input_vcf"
    exit 1
fi

if [[ ! -f "$reference" ]]; then
    echo "[ERROR] Reference FASTA not found: $reference"
    exit 1
fi

if [[ ! -f "${reference}.fai" ]]; then
    echo "[ERROR] Reference FASTA index (.fai) not found: ${reference}.fai"
    exit 1
fi

echo "[INFO] Normalizing VCF: $input_vcf"
bcftools norm -f "$reference" -m -snps -m -indels -Oz -o "$output_vcf" "$input_vcf"
tabix -p vcf "$output_vcf"
echo "[INFO] Indexed normalized VCF: $output_vcf"

#gunzip
if [[ -n "$decompressed_vcf" ]]; then
    echo "[INFO] Decompressing $output_vcf → $decompressed_vcf"
    gunzip -c "$output_vcf" > "$decompressed_vcf"
    echo "[INFO] Decompressed file written to: $decompressed_vcf"
fi

echo "[INFO] Normalization (and optional decompression) completed successfully."
