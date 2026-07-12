import pandas as pd

# ==========================================================
# Load Final Candidate Variants
# ==========================================================

df = pd.read_excel("Annotation/Final_Candidate_Variants.xlsx")

# Number of final canonical/MANE variants
canonical_variants = len(df)

# Top 2 prioritized variants (already sorted in previous script)
top_variants = df.head(2)

# ==========================================================
# Create Clinical Report
# ==========================================================

report_path = "Annotation/Clinical_Report.txt"

with open(report_path, "w") as report:

    report.write("="*54 + "\n")
    report.write("ADPKD Whole Exome Sequencing Clinical Report\n")
    report.write("="*54 + "\n\n")

    report.write("Sample:\n")
    report.write("SRR21384731\n\n")

    report.write("Reference Genome:\n")
    report.write("GRCh38\n\n")

    report.write("Pipeline:\n")
    report.write("FastQC\n")
    report.write("BWA-MEM\n")
    report.write("SAMtools\n")
    report.write("GATK\n")
    report.write("VEP\n\n")

    report.write("-"*54 + "\n")
    report.write("Summary\n")
    report.write("-"*54 + "\n\n")

    report.write("Total Variants Detected      : 154192\n")
    report.write("PKD Gene Variants            : 642\n")
    report.write("Candidate Variants           : 44\n")
    report.write(f"Canonical Candidate Variants : {canonical_variants}\n\n")

    report.write("-"*54 + "\n")
    report.write("Top Candidate Variants\n")
    report.write("-"*54 + "\n\n")

    for index, row in enumerate(top_variants.itertuples(index=False), start=1):

        report.write(f"{index}.\n\n")

        report.write("Gene:\n")
        report.write(f"{row.Symbol}\n\n")

        report.write("Variant:\n")
        report.write(f"{row.Uploaded_variation}\n\n")

        report.write("Consequence:\n")
        report.write(f"{row.Consequence}\n\n")

        report.write("Interpretation:\n")
        report.write(f"{row.Interpretation}\n\n")

        report.write("-"*54 + "\n\n")

    report.write("Figures Generated\n\n")

    report.write("Figures saved in: Figures/\n\n")

    report.write("Figure1  Interpretation Distribution\n")
    report.write("Figure2  Consequence Distribution\n")
    report.write("Figure3  Impact Distribution\n")
    report.write("Figure4  Canonical Transcript Distribution\n")
    report.write("Figure5  PKD1 vs PKD2 Distribution\n\n")

    report.write("-"*54 + "\n")
    report.write("Conclusion\n")
    report.write("-"*54 + "\n\n")

    report.write(
        "The analysis identified two high-priority candidate variants "
        "in PKD1 after functional annotation and prioritization.\n\n"
    )

    report.write(
        "One variant is a high-impact stop_gained mutation, while the "
        "other is a rare missense variant predicted to be damaging.\n\n"
    )

    report.write(
        "Additional variants are likely benign or lower-priority "
        "variants based on consequence and population frequency.\n\n"
    )

    report.write(
        "These findings should be interpreted together with ClinVar "
        "evidence, ACMG/AMP guidelines, family segregation analysis, "
        "and patient phenotype before clinical decision-making.\n\n"
    )

    report.write("="*54 + "\n")

print("="*50)
print("Clinical Report Generated Successfully!")
print("="*50)
print(f"Saved to: {report_path}")
