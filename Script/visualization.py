# ============================================================
# ADPKD Whole Exome Sequencing Pipeline
# Visualization Script
# ============================================================

import os
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Read Final Candidate Variant Table
# ============================================================

df = pd.read_excel("Annotation/Final_Candidate_Variants.xlsx")

# ============================================================
# Create Output Folder
# ============================================================

output_dir = "Figures"
os.makedirs(output_dir, exist_ok=True)

# ============================================================
# Publication Style
# ============================================================

plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300

plt.rcParams["font.size"] = 11
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10

# ============================================================
# Function to Save Bar Plots
# ============================================================

def save_barplot(series, title, xlabel, ylabel, filename):

    plt.figure(figsize=(8,6))

    ax = series.plot(kind="bar")

    plt.title(title, fontweight="bold")

    plt.xlabel(xlabel)

    plt.ylabel(ylabel)

    plt.grid(axis="y", linestyle="--", alpha=0.4)

    plt.xticks(rotation=25, ha="right")

    for container in ax.containers:
        ax.bar_label(container, padding=3)

    plt.tight_layout()

    plt.savefig(
        os.path.join(output_dir, filename),
        bbox_inches="tight"
    )

    plt.close()

# ============================================================
# Figure 1
# Interpretation Distribution
# ============================================================

interpretation_counts = df["Interpretation"].value_counts()

save_barplot(
    interpretation_counts,
    "Interpretation Distribution",
    "Interpretation",
    "Number of Variants",
    "Figure1_Interpretation_Distribution.png"
)

# ============================================================
# Figure 2
# Consequence Distribution
# ============================================================

consequence_counts = df["Consequence"].value_counts()

save_barplot(
    consequence_counts,
    "Variant Consequence Distribution",
    "Consequence",
    "Number of Variants",
    "Figure2_Consequence_Distribution.png"
)

# ============================================================
# Figure 3
# Impact Distribution
# ============================================================

impact_counts = df["Impact"].value_counts()

save_barplot(
    impact_counts,
    "Variant Impact Distribution",
    "Impact",
    "Number of Variants",
    "Figure3_Impact_Distribution.png"
)

# ============================================================
# Figure 4
# Canonical Transcript Distribution
# ============================================================

canonical = (
    df["Canonical"]
    .fillna("NO")
    .replace("", "NO")
)

canonical_counts = canonical.value_counts()

save_barplot(
    canonical_counts,
    "Canonical Transcript Distribution",
    "Canonical Transcript",
    "Number of Variants",
    "Figure4_Canonical_Transcript_Distribution.png"
)

# ============================================================
# Figure 5
# PKD1 vs PKD2 Distribution
# ============================================================

gene_counts = df["Symbol"].value_counts()

plt.figure(figsize=(6,6))

plt.pie(
    gene_counts,
    labels=gene_counts.index,
    autopct="%1.1f%%",
    startangle=90
)

plt.title(
    "PKD1 vs PKD2 Candidate Variant Distribution",
    fontweight="bold"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        output_dir,
        "Figure5_PKD1_vs_PKD2_PieChart.png"
    ),
    bbox_inches="tight"
)

plt.close()

# ============================================================
# Completed
# ============================================================

print("="*60)
print("Visualization Completed Successfully")
print("="*60)

print("\nFigures saved to:")

print(output_dir)
