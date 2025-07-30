#!/bin/bash

# Usage: bash bcftools_pass.sh <input_filtered_vcf.gz> <output_pass_vcf.gz>

set -e

# Read arguments
input_vcf="$1"
output_vcf="$2"

# Check that input exists
if [[ ! -f "$input_vcf" ]]; then
    echo "Error: input file '$input_vcf' not found!"
    exit 1
fi

# Apply bcftools filter
bcftools view -f PASS -Oz -o "$output_vcf" "$input_vcf"

# Index the output VCF
tabix -p vcf "$output_vcf"

echo "PASS-filtered VCF written to: $output_vcf"
