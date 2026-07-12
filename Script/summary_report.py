import pandas as pd

# -----------------------------
# Load Final Candidate Variants
# -----------------------------
df = pd.read_excel("Annotation/Final_Candidate_Variants.xlsx")

# -----------------------------
# Basic Statistics
# -----------------------------
total_variants = len(df)

gene_counts = df["Symbol"].value_counts()

consequence_counts = df["Consequence"].value_counts()

impact_counts = df["Impact"].value_counts()

interpretation_counts = df["Interpretation"].value_counts()

canonical_count = (df["Canonical"] == "YES").sum()

mane_count = df["MANE"].notna().sum()

# -----------------------------
# Write Report
# -----------------------------
with open("Annotation/Summary_Report.txt", "w") as report:

    report.write("=" * 60 + "\n")
    report.write("        ADPKD WES Candidate Variant Summary Report\n")
    report.write("=" * 60 + "\n\n")

    report.write(f"Total Final Candidate Variants : {total_variants}\n\n")

    report.write("Gene Distribution\n")
    report.write("-----------------\n")

    for gene, count in gene_counts.items():
        report.write(f"{gene:<10} : {count}\n")

    report.write("\n")

    report.write("Variant Consequences\n")
    report.write("--------------------\n")

    for consequence, count in consequence_counts.items():
        report.write(f"{consequence:<25} : {count}\n")

    report.write("\n")

    report.write("Impact Distribution\n")
    report.write("-------------------\n")

    for impact, count in impact_counts.items():
        report.write(f"{impact:<15} : {count}\n")

    report.write("\n")

    report.write("Interpretation Summary\n")
    report.write("----------------------\n")

    for interpretation, count in interpretation_counts.items():
        report.write(f"{interpretation:<70} : {count}\n")

    report.write("\n")

    report.write(f"Canonical Variants : {canonical_count}\n")
    report.write(f"MANE Transcripts   : {mane_count}\n")

    report.write("\n")

    report.write("=" * 60 + "\n")
    report.write("End of Report\n")
    report.write("=" * 60 + "\n")

print("Summary report generated successfully!")
