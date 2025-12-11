# ADPKD WES GATK Analysis Pipeline
## 🌡️ Pipeline Overview

This pipeline analyzes Whole Exome Sequencing (WES) data for Autosomal Dominant Polycystic Kidney Disease (ADPKD) using sample SRR21384731, focusing on chromosomes 4 and 16 that harbor the PKD1 and PKD2 genes. It performs raw read quality control, adapter and quality trimming, targeted alignment to a combined chr4/chr16 GRCh38 reference, BAM processing, BQSR, variant calling, SNP hard‑filtering and ANNOVAR‑ready output generation.

**Pipeline Workflow Summary:**

FASTQ → FastQC → FastP → BWA‑MEM → GATK (MarkDuplicates + BQSR + HaplotypeCaller) → VariantFiltration → ANNOVAR input.

# 🚀 Quick Start
## 1. Clone & Setup Environment

    git clone https://github.com/sathyaviswanath/ADPKD_WES_GATK_Pipeline.git
    cd ADPKD_WES_GATK_Pipeline

    chmod +x run_pipeline.sh

Install core tools (if not already installed) and ensure Docker is available:


    sudo apt update && sudo apt upgrade -y
    sudo apt -y install fastqc fastp bwa samtools bcftools vcftools docker.io
    sudo docker pull broadinstitute/gatk:latest
    sudo docker pull bioinfochrustrasbourg/annovar:latest

The pipeline script also attempts to install the main utilities and pull the GATK image during execution.

## 2. Download Data & Run Pipeline

All downloads and processing are handled by run_pipeline.sh:

    ./run_pipeline.sh

This will, Create Raw_Data/ and Outputs/ directories under the current folder.

- Download paired‑end FASTQ files for SRR21384731 from ENA into Raw_Data/.

- Download GRCh38 chr4 and chr16 FASTA files, build chr4_chr16.fa and run the full GATK‑based variant discovery and ANNOVAR‑prep workflow into Outputs/.

# 📁 Documentation

Detailed documentation is provided in the [Documentation/ folder](Documentation/folder):

1. See [Documentation/1.Pipeline_Overview.md](Documentation/1.Pipeline_Overview.md) for:

- Scientific context for ADPKD, description of raw data and chr4/chr16 reference, analysis goals, and high‑level workflow summary.

2. See [Documentation/2.Pipeline.md](Documentation/2.Pipeline.md) for:

- Tool specification table and Step‑by‑step pipeline description with all key commands, expected inputs/outputs at each stage.

3. See [Documentation/3.Troubleshooting.md](Documentation/3.Troubleshooting.md) for:

- Common issues (e.g., storage limits in Codespaces, Docker failures) and strategies to solve it.

# 📊 Key Outputs
All important results are written to Outputs/:

1. QC Reports
2. Alignment / BAM Files
3. Variant Files
4. Annotation Files
5. **Outputs/humandb/hg38_refGene*.gz:** hg38 refGene annotation databases used for functional annotation.

# ⚙️ Customization
To run this pipeline on a different sample or modify behavior:

- Edit the variable block near the top of run_pipeline.sh:

- Change SAMPLE_ID and SRA path if using another sample.

- Point reference URLs to different assemblies if needed.

For whole‑genome or other‑region analysis:

- Replace the chr4/chr16 FASTA URLs with full‑genome GRCh38 FASTA.

- Remove or adjust the bcftools view -r chr4,chr16 step for known sites.

# 📚 References
## Data Sources Used in This Pipeline

- **Raw WES sample (SRR21384731)**  
  
  Paired‑end FASTQ files downloaded from the European Nucleotide Archive / SRA:

  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/031/SRR21384731/SRR21384731_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/031/SRR21384731/SRR21384731_2.fastq.gz

- **GRCh38 Reference Genome (chr4 and chr16)**  

  Chromosome‑level FASTA files from NCBI GRCh38.p14 assembly:  

  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr4.fna.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr16.fna.gz

- **dbSNP Known Sites (hg38)**  
  
  dbSNP v138 VCF used for base quality score recalibration:  
  
  https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz

- **ANNOVAR hg38 refGene Databases**  
  
  Annotation databases downloaded by the pipeline:  
  
   http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
   http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
   http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz
