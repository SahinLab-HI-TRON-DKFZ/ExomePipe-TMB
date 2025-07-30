The pipeline includes the following steps:
- Variant calling using GATK Mutect2
- Orientation bias correction
- Contamination estimation
- Variant filtering
- VEP annotation
- Variant normalization
- Conversion of VCF to TSV
- TMB calculation

Run the pipeline

```
snakemake --cores 2
```
