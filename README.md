The pipeline includes the following steps:
- Variant calling using GATK Mutect2
- Orientation bias correction
- Contamination estimation
- Variant filtering
- VEP annotation
- Variant normalization
- Conversion of VCF to TSV
- TMB calculation

First, activate your Snakemake conda environment.
Run the pipeline

```
snakemake --cores 2
```
```
project/
│
├── Snakefile
├── config/
│ └── samples.tsv # Sample metadata
├── scripts/
│ ├── bcftools_pass.sh
│ ├── vep_run.sh
│ ├── bcftools_norm.sh
│ ├── vcf2tsv.sh
│ └── tmb_calc.py
├── envs/
│ └── vep.yaml 
├── reference/
└── results/
```

Author: Sakshi Singh

Email: [sakshi.singh@dkfz-heidelberg.de](mailto:sakshi.singh@dkfz-heidelberg.de)  
