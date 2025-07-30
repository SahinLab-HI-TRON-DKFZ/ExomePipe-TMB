import pandas as pd
import pysam

# Load sample metadata
samples_df = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = list(samples_df["sample"])

# Extract SM tag from BAM header
def get_sample_name(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    rg_list = bam.header.get("RG", [])
    return rg_list[0]["SM"] if rg_list else "UNKNOWN"

# Build sample info dict
SAMPLE_INFO = {
    row["sample"]: {
        "tumor_bam": row["tumor_bam"],
        "normal_bam": row["normal_bam"],
        "tumor_sample": get_sample_name(row["tumor_bam"]),
        "normal_sample": get_sample_name(row["normal_bam"])
    }
    for _, row in samples_df.iterrows()
}

# References
REF_FASTA = "reference/BWA_GRCh38/Homo_sapiens_assembly38.fasta"
GNOMAD_VCF = "reference/GnomAD_GRCh38/gnomad.joint.v4.1.NFE_AC_AF.vcf.gz"
GATK_JAR = "/home/d194-admin/master-mainz/Softwares/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar"

# Final outputs
rule all:
    input:
        expand("results/{sample}_TMBnormalized_CSQ.tsv", sample=SAMPLES),
        expand("results/{sample}_TMB_plot.png", sample=SAMPLES),
        expand("results/{sample}_TMB_summary.tsv", sample=SAMPLES)

# Mutect2 variant calling
rule mutect2:
    input:
        tumor=lambda wc: SAMPLE_INFO[wc.sample]["tumor_bam"],
        normal=lambda wc: SAMPLE_INFO[wc.sample]["normal_bam"]
    output:
        unfiltered="results/{sample}_unfiltered.vcf.gz",
        f1r2="results/{sample}_f1r2.tar.gz"
    params:
        tumor_sample=lambda wc: SAMPLE_INFO[wc.sample]["tumor_sample"],
        normal_sample=lambda wc: SAMPLE_INFO[wc.sample]["normal_sample"],
        reference=REF_FASTA,
        gnomad=GNOMAD_VCF,
        gatk=GATK_JAR
    shell:
        """
        java -jar {params.gatk} Mutect2 \
            -R {params.reference} \
            -I {input.tumor} \
            -I {input.normal} \
            --tumor-sample {params.tumor_sample} \
            --normal-sample {params.normal_sample} \
            --germline-resource {params.gnomad} \
            -O {output.unfiltered} \
            --f1r2-tar-gz {output.f1r2}
        """

# LearnReadOrientationModel
rule learn_orientation:
    input:
        f1r2="results/{sample}_f1r2.tar.gz"
    output:
        model="results/{sample}.read_orientation_model.tar.gz"
    params:
        gatk=GATK_JAR
    shell:
        """
        java -jar {params.gatk} LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.model}
        """

# GetPileupSummaries
rule pileup_summaries:
    input:
        bam=lambda wc: SAMPLE_INFO[wc.sample]["tumor_bam"],
        gnomad=GNOMAD_VCF
    output:
        table="results/{sample}.pileup_summary.table"
    params:
        gatk=GATK_JAR
    shell:
        """
        java -jar {params.gatk} GetPileupSummaries \
            -I {input.bam} \
            --variant {input.gnomad} \
            --intervals {input.gnomad} \
            --minimum-population-allele-frequency 0.05 \
            -O {output.table}
        """

# Dummy tumor segmentation table (replace with real segmentation rule if needed)
rule dummy_segments:
    output:
        "results/{sample}.segments.table"
    shell:
        "touch {output}"

# CalculateContamination
rule calculate_contamination:
    input:
        pileup="results/{sample}.pileup_summary.table",
        segments="results/{sample}.segments.table"
    output:
        contamination="results/{sample}.contamination.table"
    params:
        gatk=GATK_JAR
    shell:
        """
        java -jar {params.gatk} CalculateContamination \
            -I {input.pileup} \
            --tumor-segmentation {input.segments} \
            -O {output.contamination}
        """

# FilterMutectCalls
rule filter_calls:
    input:
        vcf="results/{sample}_unfiltered.vcf.gz",
        contamination="results/{sample}.contamination.table",
        segments="results/{sample}.segments.table",
        model="results/{sample}.read_orientation_model.tar.gz"
    output:
        filtered="results/{sample}_mutect2.filtered.vcf.gz"
    params:
        gatk=GATK_JAR,
        reference=REF_FASTA
    shell:
        """
        java -jar {params.gatk} FilterMutectCalls \
            -R {params.reference} \
            -V {input.vcf} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --ob-priors {input.model} \
            -O {output.filtered}
        """

# PASS filtering
rule filter_pass:
    input:
        "results/{sample}_mutect2.filtered.vcf.gz"
    output:
        "results/{sample}.mutect2.filtered.PASS.vcf.gz"
    shell:
        "bash scripts/bcftools_pass.sh {input} {output}"

# VEP annotation
rule vep:
    input:
        "results/{sample}.mutect2.filtered.PASS.vcf.gz"
    output:
        "results/{sample}_VEP_TMBfiltered.vcf.gz"
    params:
        sample="{sample}",
        outdir="results"
    conda:
        "envs/vep.yaml"
    shell:
        "bash scripts/vep_run.sh {input} {params.sample} {params.outdir}"

# Normalization
rule normalize:
    input:
        vcf="results/{sample}_VEP_TMBfiltered.vcf.gz",
        ref=REF_FASTA
    output:
        norm="results/{sample}_normalized.vcf.gz",
        decompressed="results/{sample}_normalized.vcf"
    shell:
        "bash scripts/bcftools_norm.sh {input.vcf} {output.norm} {input.ref} {output.decompressed}"

# Convert VCF to TSV
rule vcf2tsv:
    input:
        "results/{sample}_normalized.vcf"
    output:
        "results/{sample}_TMBnormalized.tsv"
    shell:
        "bash scripts/vcf2tsv.sh {input} {output}"

# TMB calculation
rule tmb_calc:
    input:
        "results/{sample}_TMBnormalized.tsv"
    output:
        csq="results/{sample}_TMBnormalized_CSQ.tsv",
        plot="results/{sample}_TMB_plot.png",
        summary="results/{sample}_TMB_summary.tsv"
    shell:
        "python scripts/tmb_calc.py {input} {output.csq} {output.plot} {output.summary}"
