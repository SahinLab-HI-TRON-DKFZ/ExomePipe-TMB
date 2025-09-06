#!/usr/bin/env python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

# Usage: python tmb_calc.py <input_tsv> <output_csq_tsv> <output_plot_png> <summary_tsv>

if len(sys.argv) != 5:
    print("Usage: python tmb_calc.py <input_TMBnormalized.tsv> <output_CSQ.tsv> <output_plot.png> <summary_TMB.tsv>")
    sys.exit(1)

input_tsv = sys.argv[1]
output_csq = sys.argv[2]
output_plot = sys.argv[3]
summary_tsv = sys.argv[4]

EXOME_SIZE_MB = 60  # Exome size (SureSelect Human All Exon V6)

df = pd.read_csv(input_tsv, sep="\t", dtype=str).fillna('')

def split_annotation_column(df, column_name, prefix):
    if column_name in df.columns:
        df[column_name] = df[column_name].astype(str).fillna('')
        expanded = df[column_name].str.split(",", expand=True)
        expanded = expanded.applymap(lambda x: x.split("|") if isinstance(x, str) else [""] * 20)
        max_cols = max(expanded.applymap(len).max())
        expanded = expanded.apply(lambda row: pd.Series(row[0] if isinstance(row[0], list) else [""] * max_cols), axis=1)
        expanded.columns = [f"{prefix}_{i}" for i in range(len(expanded.columns))]
        df = df.drop(columns=[column_name]).join(expanded)
    return df

df = split_annotation_column(df, "CSQ", "CSQ")
df.to_csv(output_csq, sep="\t", index=False)
print(f"[INFO] CSQ-processed TSV saved: {output_csq}")

# === TMB Calculation ===
df = pd.read_csv(output_csq, sep="\t")
if "CSQ_1" not in df.columns:
    print(f"[ERROR] CSQ_1 column not found in {output_csq}")
    sys.exit(1)

mutation_counts = df["CSQ_1"].value_counts()
colors = [to_hex(c) for c in plt.cm.tab20.colors[:len(mutation_counts)]]

total_mutations = mutation_counts.sum()
tmb_exome = total_mutations / EXOME_SIZE_MB
fig, ax = plt.subplots(figsize=(12, 8))
mutation_counts.plot(kind="bar", color=colors, edgecolor='black', linewidth=0.8, ax=ax)

ax.set_title(f"TMB Score", fontsize=10)
ax.set_xlabel("Annotation Type", fontsize=12)
ax.set_ylabel("Counts", fontsize=12)
ax.set_xticklabels(mutation_counts.index, rotation=45, ha="right")

for i, v in enumerate(mutation_counts):
    ax.text(i, v + max(mutation_counts) * 0.02, str(v), ha="center", va="bottom", fontsize=10, fontweight="bold")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

tmb_text = f"Total Mutations: {total_mutations}\nTMB (Exome {EXOME_SIZE_MB} Mb): {tmb_exome:.2f} mutations/Mb"
props = dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="black", alpha=0.9)
ax.text(len(mutation_counts) - 0.5, max(mutation_counts) * 0.8, tmb_text, fontsize=12, bbox=props, ha="left", va="top")

plt.tight_layout()
plt.savefig(output_plot)
print(f"[INFO] TMB plot saved: {output_plot}")

df_summary = pd.DataFrame([{
    "Sample Name": os.path.basename(input_tsv).replace("_TMBnormalized.tsv", ""),
    "Total Mutations": total_mutations,
    "TMB Score": round(tmb_exome, 2)
}])
df_summary.to_csv(summary_tsv, sep="\t", index=False)
print(f"[INFO] TMB summary written: {summary_tsv}")
