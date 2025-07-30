#!/bin/bash

# Usage:
# bash vep_run.sh <input_vcf.gz> <sample_id> <output_dir>

set -e

# === Input Arguments ===
INPUT_VCF="$1"
SAMPLE_ID="$2"
OUTDIR="$3"

# === Configuration ===
REFERENCE_FASTA="reference/BWAIndex/Homo_sapiens_assembly38.fasta"
VEP_CACHE_DIR="/home/d194-admin/vep_data"
VEP_CACHE_VERSION="113"
ASSEMBLY="GRCh38"

# === Output Files ===
VEP_OUT_VCF="${OUTDIR}/${SAMPLE_ID}_VEP.vcf.gz"
VEP_HTML="${OUTDIR}/${SAMPLE_ID}_VEP_stats.html"
VEP_SYN_VCF="${OUTDIR}/${SAMPLE_ID}_VEP_removesynonymous.vcf.gz"
VEP_TMB_VCF="${OUTDIR}/${SAMPLE_ID}_VEP_TMBfiltered.vcf.gz"

# === Validation Checks ===
[[ ! -f "$REFERENCE_FASTA" ]] && echo "[ERROR] Missing reference FASTA: $REFERENCE_FASTA" && exit 1
[[ ! -f "${REFERENCE_FASTA}.fai" ]] && echo "[ERROR] Missing FASTA index: ${REFERENCE_FASTA}.fai" && exit 1
[[ ! -d "${VEP_CACHE_DIR}/homo_sapiens/${VEP_CACHE_VERSION}_${ASSEMBLY}" ]] && \
    echo "[ERROR] VEP cache not found at ${VEP_CACHE_DIR}/homo_sapiens/${VEP_CACHE_VERSION}_${ASSEMBLY}" && exit 1

mkdir -p "$OUTDIR"

# === Run VEP Annotation ===
echo "[INFO] Running VEP annotation on $INPUT_VCF..."

vep --assembly "$ASSEMBLY" \
    --cache \
    --dir_cache "$VEP_CACHE_DIR" \
    --cache_version "$VEP_CACHE_VERSION" \
    --fasta "$REFERENCE_FASTA" \
    --fork 6 \
    --format vcf \
    --compress_output bgzip \
    --vcf \
    --hgvs \
    --symbol \
    --tsl \
    --appris \
    --canonical \
    --check_existing \
    --pick \
    --protein \
    --check_ref \
    --offline \
    -i "$INPUT_VCF" \
    -o "$VEP_OUT_VCF" \
    --stats_file "$VEP_HTML" || { echo "[ERROR] VEP failed"; exit 1; }

# === Filtering Step ===
echo "[INFO] Removing synonymous variants..."
bcftools filter -e 'INFO/CSQ ~ "synonymous_variant"' \
    "$VEP_OUT_VCF" -Oz -o "$VEP_SYN_VCF"
tabix -p vcf "$VEP_SYN_VCF"

echo "[INFO] Filtering for TMB-relevant variants..."
bcftools filter -e 'INFO/CSQ ~ "synonymous_variant" || \
INFO/CSQ ~ "5_prime_UTR_variant" || INFO/CSQ ~ "3_prime_UTR_variant" || \
INFO/CSQ ~ "intron_variant" || INFO/CSQ ~ "upstream_gene_variant" || \
INFO/CSQ ~ "downstream_gene_variant" || INFO/CSQ ~ "TFBS_ablation" || \
INFO/CSQ ~ "TFBS_amplification" || INFO/CSQ ~ "TF_binding_site_variant" || \
INFO/CSQ ~ "regulatory_region_ablation" || INFO/CSQ ~ "regulatory_region_amplification" || \
INFO/CSQ ~ "regulatory_region_variant" || INFO/CSQ ~ "intergenic_variant"' \
    "$VEP_OUT_VCF" -Oz -o "$VEP_TMB_VCF"
tabix -p vcf "$VEP_TMB_VCF"

# === Summary ===
echo "[INFO] VEP annotation and filtering completed for $SAMPLE_ID"
echo "Annotated VCF:         $VEP_OUT_VCF"
echo "No-synonymous VCF:     $VEP_SYN_VCF"
echo "TMB-filtered VCF:      $VEP_TMB_VCF"
echo "VEP HTML stats report: $VEP_HTML"
