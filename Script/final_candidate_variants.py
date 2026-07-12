import pandas as pd

# Read VEP Annotation file

vep = pd.read_csv(
	"Annotation/SRR21384731_filtered_snps_annotated.txt",
	sep="\t",
	comment="#"
)

# Display basic information
print("Shape: ", vep.shape)

# Assign Column Names
vep.columns = [
	"Uploaded_variation",
	"Location",
	"Allele",
	"Gene",
	"Feature",
	"Feature_type",
	"Consequence",
	"cDNA_position",
	"CDS_position",
	"Protein_position",
	"Amino_acids",
	"Codons",
	"Existing_variation",
	"Extra"
]
print(vep.head())

# Function to extract values from Extra column
def extract_field(extra, field):
	for item in str(extra).split(";"):
		if item.startswith(field + "="):
			return item.split("=", 1)[1]
	return None

vep["Impact"] = vep["Extra"].apply(lambda x: extract_field(x, "IMPACT"))
vep["Symbol"] = vep["Extra"].apply(lambda x: extract_field(x, "SYMBOL"))
vep["ClinVar"] = vep["Extra"].apply(lambda x: extract_field(x, "CLIN_SIG"))
vep["SIFT"] = vep["Extra"].apply(lambda x: extract_field(x, "SIFT"))
vep["PolyPhen"] = vep["Extra"].apply(lambda x: extract_field(x, "PolyPhen"))
vep["gnomAD_AF"] = vep["Extra"].apply(lambda x: extract_field(x, "MAX_AF"))
vep["Canonical"] = vep["Extra"].apply(lambda x: extract_field(x, "CANONICAL"))
vep["MANE"] = vep["Extra"].apply(lambda x: extract_field(x, "MANE_SELECT"))
print(vep[[
	"Uploaded_variation",
	"Consequence",
        "Impact",
        "Symbol",
        "ClinVar",
        "SIFT",
        "PolyPhen",
        "gnomAD_AF",
        "Canonical",
        "MANE"
	]].head()
)

# Keep only PKD1 and PKD2 variants
pkd = vep[vep["Symbol"].isin(["PKD1", "PKD2"])]
print("\nTotal PKD Variants:", len(pkd))
print(pkd[["Uploaded_variation", "Symbol", "Impact"]].head())

# Keep HIGH and MODERATE impact variants
candidate = pkd[pkd["Impact"].isin(["HIGH", "MODERATE"])]
print("\nCandidate Variants:", len(candidate))
print(candidate[["Uploaded_variation", "Symbol", "Impact", "Consequence"]].head(10))

# Keep canonical transcripts if available
candidate = candidate[(candidate["Canonical"] == "YES") | (candidate["MANE"].notna())]
print("\nCanonical/MANE candidate variants:", len(candidate))

# Select important columns for the final table
final_table = candidate[
    [
        "Uploaded_variation",
        "Location",
        "Symbol",
        "Consequence",
        "Impact",
        "ClinVar",
        "SIFT",
        "PolyPhen",
        "gnomAD_AF",
        "Canonical",
        "MANE"
    ]
]
print("\nFinal Candidate Variants:")
print(final_table)

# Function to interpret candidate variants
def interpret_variant(row):

    consequence = str(row["Consequence"])
    polyphen = str(row["PolyPhen"])
    af = row["gnomAD_AF"]

    # Convert gnomAD_AF to float if possible
    try:
        af = float(af)
    except:
        af = None

    # High-impact loss-of-function variants
    if "stop_gained" in consequence:
        return "High-impact loss-of-function variant; Strong Candidate Variant"

    # Missense variants
    elif "missense_variant" in consequence:

        if af is not None and af < 0.01:

            if "probably_damaging" in polyphen:
                return "Rare missense variant predicted damaging; High-priority candidate"

            elif "possibly_damaging" in polyphen:
                return "Rare missense variant possibly damaging; Requires further evaluation"

            elif "benign" in polyphen:
                return "Rare missense variant but predicted benign"

            else:
                return "Rare missense variant; Functional significance uncertain"

        elif af is not None and af >= 0.01:
            return "Common population variant; Low priority for further investigation"

        else:
            return "Missense variant; Insufficient population frequency data"

    # Other variants
    else:
        return "Lower-priority variant"


# Add interpretation column
final_table["Interpretation"] = final_table.apply(interpret_variant, axis=1)
print(final_table[[
    "Uploaded_variation",
    "Location",
    "Symbol",
    "Consequence",
    "Impact",
    "PolyPhen",
    "gnomAD_AF",
    "Canonical",
    "MANE",
    "Interpretation"
]])
print("\nInterpretation added to final_table successfully!!")

# Sorting interpretation column from Strong candidate variant to lower-priority variant
priority = {
    "High-impact loss-of-function variant; Strong Candidate Variant": 1,
    "Rare missense variant predicted damaging; High-priority candidate": 2,
    "Rare missense variant possibly damaging; Requires further information": 3,
    "Rare missense variant but predicted benign": 4,
    "Rare missense variant; Functional significance uncertain": 5,
    "Missense variant; insufficient population frequency data": 6,
    "Common population missense variant; less likely disease-causing": 7,
    "Lower-priority variant": 8
}

final_table["Priority"] = final_table["Interpretation"].map(priority)
final_table = final_table.sort_values(by=["Priority", "Symbol", "Location"])

# Removing priority column
final_table = final_table.drop(columns="Priority")

print("\nSorted Final Candidate Variants:")
print(final_table)

# Save to excel
final_table.to_excel("Annotation/Final_Candidate_Variants.xlsx", index=False)
print("\nFinal candidate variant table saved successfully!")
