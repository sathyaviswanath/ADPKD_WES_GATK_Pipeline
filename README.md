# 🧬 ADPKD Whole Exome Sequencing (WES) Analysis Pipeline

## 📖 Overview

This repository presents a reproducible **Whole Exome Sequencing (WES)** analysis pipeline developed to identify and prioritize clinically relevant variants associated with **Autosomal Dominant Polycystic Kidney Disease (ADPKD)**.

The workflow performs quality control, read alignment, variant calling, functional annotation using **Ensembl VEP**, candidate variant prioritization, visualization, summary statistics generation, and automated clinical report generation.

---

# 🩺 Autosomal Dominant Polycystic Kidney Disease (ADPKD)

Autosomal Dominant Polycystic Kidney Disease (ADPKD) is the **most common inherited kidney disease**, characterized by the progressive formation of numerous fluid-filled cysts in the kidneys and, occasionally, other organs. As cysts enlarge over time, they can lead to kidney enlargement, chronic pain, hypertension, reduced kidney function, and ultimately **kidney failure**, affecting nearly **50% of patients by the age of 60 years**.

---

# 🧬 Genetic Basis

ADPKD is primarily caused by pathogenic variants in two genes:

- **PKD1** (Chromosome 16) – responsible for approximately **85%** of cases.
- **PKD2** (Chromosome 4) – responsible for most of the remaining cases.

These genes encode the proteins **Polycystin-1 (PC1)** and **Polycystin-2 (PC2)**, which regulate tubular epithelial cell function, calcium signaling, and primary cilia-mediated signaling pathways.

Mutations in these genes disrupt normal cellular signaling, resulting in abnormal cell proliferation, cyst formation, and the gradual decline of kidney function.

---

# 🎯 Project Objectives

- Perform quality assessment of WES data.
- Align sequencing reads to the **GRCh38** reference genome.
- Identify high-confidence germline variants using **GATK Best Practices**.
- Functionally annotate variants using **Ensembl VEP**.
- Prioritize clinically relevant variants in **PKD1** and **PKD2**.
- Generate publication-ready visualizations.
- Produce automated summary and clinical reports.

---

# 📂 Project Information

**Sample ID**

SRR21384731

**Reference Genome**

GRCh38

**Genes of Interest**

- PKD1
- PKD2

---

# ⚙️ Pipeline Workflow

```
FASTQ
   │
FastQC
   │
FastP
   │
BWA-MEM
   │
SAMtools
   │
GATK
   │
Analysis-ready SNP VCF
   │
Ensembl VEP Annotation
   │
Python Variant Prioritization
   │
Summary Report
   │
Visualization
   │
Clinical Report
```

---

# 🚀 Features

- ✅ FastQC quality assessment
- ✅ Adapter trimming using FastP
- ✅ Read alignment using BWA-MEM
- ✅ BAM processing using SAMtools
- ✅ Variant calling using GATK HaplotypeCaller
- ✅ Functional annotation using Ensembl VEP
- ✅ PKD1 and PKD2 candidate variant prioritization
- ✅ Canonical and MANE transcript selection
- ✅ Automated variant interpretation
- ✅ Summary statistics generation
- ✅ Publication-quality visualizations
- ✅ Automated clinical report generation

---

# 🔬 Downstream Analysis Strategy

Although **GATK** identifies both **Single Nucleotide Polymorphisms (SNPs)** and **small insertions/deletions (Indels)**, this project focuses on **SNP prioritization** during downstream analysis.

### Why only SNPs?

- SNPs constitute the majority of high-confidence germline variants detected in Whole Exome Sequencing.
- Most well-characterized pathogenic ADPKD variants reported in **PKD1** and **PKD2** are SNPs.
- Indels often require additional validation because repetitive genomic regions, particularly within **PKD1**, can introduce alignment and variant-calling artifacts.
- Restricting downstream analysis to SNPs provides a robust, reproducible, and clinically interpretable demonstration of the variant prioritization workflow.

Future versions of this pipeline can be extended to include comprehensive **Indel prioritization**.

---

# 📊 Final Results

| Metric | Count |
|---------|------:|
| Total Variants Detected | **154,192** |
| PKD Gene Variants | **642** |
| Candidate Variants | **44** |
| Canonical Candidate Variants | **7** |
| High-Priority Candidate Variants | **2** |

---

# 📁 Repository Structure

```
ADPKD_WES_GATK_Pipeline/

├── Annotation/
│   ├── Final_Candidate_Variants.xlsx
│   ├── Summary_Report.txt
│   ├── Clinical_Report.txt
│   └── ...
│
├── Documentation/
│   ├── 1.Pipeline_Overview.md
│   ├── 2.Pipeline.md
│   └── 3.Troubleshooting.md
│
├── Figures/
│   ├── Figure1_Interpretation_Distribution.png
│   ├── Figure2_Consequence_Distribution.png
│   ├── Figure3_Impact_Distribution.png
│   ├── Figure4_Canonical_Transcript_Distribution.png
│   └── Figure5_PKD1_vs_PKD2_PieChart.png
│
├── Outputs/
│
├── Script/
│   ├── run_pipeline.sh
│   ├── final_candidate_variants.py
│   ├── summary_report.py
│   ├── visualization.py
│   └── clinical_report.py
│
├── requirements.txt
└── README.md
```

---

---

# ▶️ Running the Pipeline

Clone the repository:

```bash
git clone https://github.com/sathyaviswanath/ADPKD_WES_GATK_Pipeline.git

cd ADPKD_WES_GATK_Pipeline
```

Make the pipeline executable:

```bash
chmod +x Script/run_pipeline.sh
```

Run the complete workflow:

```bash
bash Script/run_pipeline.sh
```

---

# 📄 Generated Outputs

The pipeline automatically generates:

- 📄 Analysis-ready SNP VCF
- 📄 VEP annotated variant file
- 📊 Final candidate variant table
- 📈 Summary statistics report
- 🏥 Clinical report
- 📉 Publication-quality figures

---

# 🔮 Future Improvements

- ACMG/AMP automated variant classification
- Comprehensive Indel prioritization
- ClinVar pathogenicity integration
- Multi-sample cohort analysis
- Interactive HTML reports
- Workflow automation using Snakemake or Nextflow

---

# 🙏 Acknowledgements

This project utilizes the following open-source tools:

- GATK
- Ensembl VEP
- FastQC
- FastP
- BWA-MEM
- SAMtools
- Python (Pandas, Matplotlib, OpenPyXL)

---

# 📜 Author

**Sathya**
*Bioinformatics Analyst*
