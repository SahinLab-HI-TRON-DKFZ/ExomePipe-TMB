#!/bin/bash

# Usage: bash mutect2_run.sh <tumor_bam> <normal_bam> <sample_id> <output_dir>

set -e
tumor_bam=$1
normal_bam=$2
sample_id=$3
outdir=$4

# Define variables
JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
GATK="/home/d194-admin/master-mainz/Softwares/gatk-4.6.0.0/gatk"
REFERENCE="reference/hs37d5_PhiX.fa"
GERMLINE="reference/gnomad.exomes.r2.1.1.sites.vcf.bgz"

mkdir -p $outdir

$GATK --java-options "-Xmx32g -XX:ParallelGCThreads=10" Mutect2 \
    -R $REFERENCE -I $tumor_bam -I $normal_bam \
    --normal-sample sample_${sample_id}_normal \
    --tumor-sample sample_${sample_id}_tumor \
    --germline-resource $GERMLINE \
    -O $outdir/${sample_id}_mutect2_unfiltered.vcf.gz \
    --f1r2-tar-gz $outdir/${sample_id}_f1r2.tar.gz

$GATK --java-options "-Xmx32g" LearnReadOrientationModel \
    -I $outdir/${sample_id}_f1r2.tar.gz \
    -O $outdir/${sample_id}_orientation_model.tar.gz

$GATK --java-options "-Xmx32g" GetPileupSummaries \
    -I $tumor_bam --variant $GERMLINE --intervals $GERMLINE \
    -O $outdir/${sample_id}_pileup_summary.table

$GATK --java-options "-Xmx32g" CalculateContamination \
    -I $outdir/${sample_id}_pileup_summary.table \
    -O $outdir/${sample_id}_contamination.table \
    --tumor-segmentation $outdir/${sample_id}_segments.table

$GATK --java-options "-Xmx32g" FilterMutectCalls \
    -R $REFERENCE \
    -V $outdir/${sample_id}_mutect2_unfiltered.vcf.gz \
    --contamination-table $outdir/${sample_id}_contamination.table \
    --tumor-segmentation $outdir/${sample_id}_segments.table \
    --ob-priors $outdir/${sample_id}_orientation_model.tar.gz \
    -O $outdir/${sample_id}_mutect2.filtered.vcf.gz
