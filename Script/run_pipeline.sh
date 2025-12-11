#!/usr/bin/env bash
set -euo pipefail

echo "ADPKD WES GATK Pipeline..."

# User-configurable variables
# Sample and read IDs
SAMPLE_ID="SRR21384731"
READ1_ID="${SAMPLE_ID}_1.fastq.gz"
READ2_ID="${SAMPLE_ID}_2.fastq.gz"

# Directories
RAW_DIR="Raw_Data"
OUT_DIR="Outputs"
HUMANDB_DIR="${OUT_DIR}/humandb"

# Reference (GRCh38)
REF_BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA"
CHR4_FNA="chr4.fna.gz"
CHR16_FNA="chr16.fna.gz"

# SRA base URL
SRA_BASE_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/031/${SAMPLE_ID}"

# Combined reference
COMBINED_REF_RAW="chr4_chr16.fa"
COMBINED_REF="${OUT_DIR}/chr4_chr16.fa"

# Known sites (dbSNP)
DBSNP_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
DBSNP_VCF="${OUT_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
KNOWN_SITES_RAW="${OUT_DIR}/known_sites_chr4_chr16.vcf.gz"
KNOWN_SITES_NC="${OUT_DIR}/known_sites_chr4_chr16_NC.vcf.gz"
CONTIG_RENAME="${OUT_DIR}/contig_rename.txt"

# Docker images
GATK_IMG="broadinstitute/gatk:latest"
ANNOVAR_IMG="bioinfochrustrasbourg/annovar:latest"

# Trimming outputs
TRIM_R1="${OUT_DIR}/R1_trimmed.fastq.gz"
TRIM_R2="${OUT_DIR}/R2_trimmed.fastq.gz"
FASTP_JSON="${OUT_DIR}/fastp.json"
FASTP_HTML="${OUT_DIR}/fastp.html"

# Alignment outputs
SAM_FILE="${OUT_DIR}/${SAMPLE_ID}.paired.sam"
DEDUP_BAM="${OUT_DIR}/${SAMPLE_ID}_sorted_dedup_reads.bam"
RECAL_TABLE="${OUT_DIR}/${SAMPLE_ID}_recal_data.table"
BQSR_BAM="${OUT_DIR}/${SAMPLE_ID}_sorted_dedup_BQSR.bam"

# Variant files
RAW_VCF="${OUT_DIR}/${SAMPLE_ID}_raw_variants.vcf"
RAW_SNPS_VCF="${OUT_DIR}/${SAMPLE_ID}_raw_snps.vcf"
RAW_INDELS_VCF="${OUT_DIR}/${SAMPLE_ID}_raw_indels.vcf"
FILTERED_SNPS_VCF="${OUT_DIR}/${SAMPLE_ID}_filtered_snps.vcf"
ANALYSIS_SNPS_VCF="${OUT_DIR}/${SAMPLE_ID}_analysis_ready_snps.vcf"
ANNOVAR_INPUT="${OUT_DIR}/${SAMPLE_ID}_filtered_snps.avinput"


# 1. Environment setup
echo "Settingup environment..."
sudo apt update
sudo apt upgrade -y
sudo apt -y install fastqc fastp bwa samtools bcftools vcftools
sudo docker pull "${GATK_IMG}"

# 2. Directories & data download
echo "Creating directories and Data Download..."
mkdir -p "${RAW_DIR}" "${OUT_DIR}"
cd "${RAW_DIR}"
wget -nc "${SRA_BASE_URL}/${READ1_ID}"
wget -nc "${SRA_BASE_URL}/${READ2_ID}"
wget -nc "${REF_BASE_URL}/${CHR4_FNA}"
wget -nc "${REF_BASE_URL}/${CHR16_FNA}"

echo "Decompress and combining reference..."
gunzip -f "./${CHR4_FNA}" "./${CHR16_FNA}"
cat chr4.fna chr16.fna > "${COMBINED_REF_RAW}"

# 3. QC on raw reads
echo "Running FASTQC..."
fastqc -o "../${OUT_DIR}" "./${READ1_ID}" "./${READ2_ID}"

# 4. Trimming and post-trim QC
echo "Running FASTP & Post FASTQC..."
fastp \
  -i "./${READ1_ID}" \
  -I "./${READ2_ID}" \
  -o "${TRIM_R1}" \
  -O "${TRIM_R2}" \
  --detect_adapter_for_pe \
  -j "${FASTP_JSON}" \
  -h "${FASTP_HTML}"

fastqc -o "../${OUT_DIR}" "${TRIM_R1}" "${TRIM_R2}"

# 5. Reference copy & indexing
echo "Indexing reference genome..."
cp "${COMBINED_REF_RAW}" "../${COMBINED_REF}"
cd "../${OUT_DIR}"

bwa index "${COMBINED_REF}"

echo "Creating Sequence Dictionary with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk CreateSequenceDictionary \
  -R /data/$(basename "${COMBINED_REF}") \
  -O /data/chr4_chr16.dict

# 6. Alignment
echo "Align reads using BWA MEM..."
bwa mem -t 4 \
  -R "@RG\tID:${SAMPLE_ID}\tPL:ILLUMINA\tSM:${SAMPLE_ID}" \
  "${COMBINED_REF}" "${TRIM_R1}" "${TRIM_R2}" \
  > "${SAM_FILE}"

# 7. Mark duplicates
echo "Marking duplicates with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk MarkDuplicatesSpark \
  -I /data/$(basename "${SAM_FILE}") \
  -O /data/$(basename "${DEDUP_BAM}") 

# 8. Known-sites VCF (dbSNP subset)
echo "Downloading & Indexing Known sites..."
wget -O "${DBSNP_VCF}" "${DBSNP_URL}"
bcftools index "${DBSNP_VCF}"

bcftools view -r chr4,chr16 "${DBSNP_VCF}" -Oz -o "${KNOWN_SITES_RAW}"
bcftools index "${KNOWN_SITES_RAW}"

# Create contig rename file
echo "Renaming VCF chromosome labels & indexing known label sites..."
cat > "${CONTIG_RENAME}" <<EOF
chr4    NC_000004.12
chr16   NC_000016.10
EOF

bcftools annotate \
  --rename-chrs "${CONTIG_RENAME}" \
  -Oz -o "${KNOWN_SITES_NC}" \
  "${KNOWN_SITES_RAW}"

tabix -p vcf "${KNOWN_SITES_NC}"

# 9. Base recalibration (BQSR)
echo "Running BaseRecalibrator with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk BaseRecalibrator \
  -I /data/$(basename "${DEDUP_BAM}") \
  -R /data/$(basename "${COMBINED_REF}") \
  --known-sites /data/$(basename "${KNOWN_SITES_NC}") \
  -O /data/$(basename "${RECAL_TABLE}")

echo "Applying BQSR with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk ApplyBQSR \
  -I /data/$(basename "${DEDUP_BAM}") \
  -R /data/$(basename "${COMBINED_REF}") \
  --bqsr-recal-file /data/$(basename "${RECAL_TABLE}") \
  -O /data/$(basename "${BQSR_BAM}")

# 10. Variant calling
echo "Calling variants with GATK HaplotypeCaller..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk HaplotypeCaller \
  -R /data/$(basename "${COMBINED_REF}") \
  -I /data/$(basename "${BQSR_BAM}") \
  -O /data/$(basename "${RAW_VCF}")

# 11. Split SNPs and INDELs
echo "Filtering SNPs and INDELs with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk SelectVariants \
  -R /data/$(basename "${COMBINED_REF}") \
  -V /data/$(basename "${RAW_VCF}") \
  --select-type SNP \
  -O /data/$(basename "${RAW_SNPS_VCF}")

sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk SelectVariants \
  -R /data/$(basename "${COMBINED_REF}") \
  -V /data/$(basename "${RAW_VCF}") \
  --select-type INDEL \
  -O /data/$(basename "${RAW_INDELS_VCF}")

# 12. Hard filter SNPs
echo "Applying variant filters with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk VariantFiltration \
  -R /data/$(basename "${COMBINED_REF}") \
  -V /data/$(basename "${RAW_SNPS_VCF}") \
  -O /data/$(basename "${FILTERED_SNPS_VCF}") \
  -filter-name "QD_filter"              -filter "QD < 2.0" \
  -filter-name "FS_filter"              -filter "FS > 60.0" \
  -filter-name "MQ_filter"              -filter "MQ < 40.0" \
  -filter-name "SOR_filter"             -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter"       -filter "MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter"  -filter "ReadPosRankSum < -8.0" \
  -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
  -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

# 13. Keep only PASS SNPs
echo "Selecting variants that pass filters with GATK..."
sudo docker run --rm -v "$PWD":/data "${GATK_IMG}" \
  gatk SelectVariants \
  --exclude-filtered \
  -V /data/$(basename "${FILTERED_SNPS_VCF}") \
  -O /data/$(basename "${ANALYSIS_SNPS_VCF}")

# 14. Convert to ANNOVAR input
echo "Converting VCF to ANNOVAR format..."
sudo docker run -it --rm -v "$PWD":/data "${ANNOVAR_IMG}" \
  perl ./convert2annovar.pl \
  -format vcf4 /data/$(basename "${FILTERED_SNPS_VCF}") \
  -outfile /data/$(basename "${ANNOVAR_INPUT}")

# 15. Download ANNOVAR databases
echo "Downloading ANNOVAR reference files..."
mkdir -p "${HUMANDB_DIR}"
cd "${HUMANDB_DIR}"

wget http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
wget http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz

echo "ADPKD WES GATK Pipeline Completed Successfully!"