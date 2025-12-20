#!/usr/bin/env bash

set -euo pipefail

echo "ADPKD WES GATK Pipeline..."

echo "Setting up environment..."
sudo apt update
sudo apt upgrade -y
sudo apt -y install fastqc fastp bwa samtools bcftools vcftools
sudo docker pull "broadinstitute/gatk:latest"

echo "Creating directories and Data Download..."
mkdir -p "Raw_Data" "Outputs"
cd "Raw_Data"
wget -nc "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/031/SRR21384731/SRR21384731_1.fastq.gz"
wget -nc "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/031/SRR21384731/SRR21384731_2.fastq.gz"
wget -nc "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr4.fna.gz"
wget -nc "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr16.fna.gz"

echo "Decompress and combining reference..."
gunzip -f "./chr4.fna.gz" "./chr16.fna.gz"
cat chr4.fna chr16.fna > "Outputs/chr4_chr16.fa"

echo "Running FASTQC..."
fastqc -o "../Outputs" "./SRR21384731_1.fastq.gz" "./SRR21384731_2.fastq.gz"

echo "Running FASTP & Post FASTQC..."
fastp \
  -i "./SRR21384731_1.fastq.gz" \
  -I "./SRR21384731_2.fastq.gz" \
  -o "Outputs/R1_trimmed.fastq.gz" \
  -O "Outputs/R2_trimmed.fastq.gz" \
  --detect_adapter_for_pe \
  -j "Outputs/fastp.json" \
  -h "Outputs/fastp.html"
fastqc -o "../Outputs" "Outputs/R1_trimmed.fastq.gz" "Outputs/R2_trimmed.fastq.gz"

echo "Indexing reference genome..."
cd "../Outputs"
bwa index "chr4_chr16.fa"

echo "Creating Sequence Dictionary with GATK..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk CreateSequenceDictionary \
  -R /data/chr4_chr16.fa \
  -O /data/chr4_chr16.dict

echo "Align reads using BWA MEM..."
bwa mem -t 4 \
  -R "@RG\tID:SRR21384731\tPL:ILLUMINA\tSM:SRR21384731" \
  "chr4_chr16.fa" "R1_trimmed.fastq.gz" "R2_trimmed.fastq.gz" \
  > "SRR21384731.paired.sam"

echo "Marking duplicates with GATK..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk MarkDuplicatesSpark \
  -I /data/SRR21384731.paired.sam \
  -O /data/SRR21384731_sorted_dedup_reads.bam

echo "Downloading & Indexing Known sites..."
wget -O "Homo_sapiens_assembly38.dbsnp138.vcf.gz" "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
bcftools index "Homo_sapiens_assembly38.dbsnp138.vcf.gz"
bcftools view -r chr4,chr16 "Homo_sapiens_assembly38.dbsnp138.vcf.gz" -Oz -o "known_sites_chr4_chr16.vcf.gz"
bcftools index "known_sites_chr4_chr16.vcf.gz"

echo "Renaming VCF chromosome labels & indexing known label sites..."
cat > "contig_rename.txt" << EOF
chr4	chr4
chr16	chr16
EOF
bcftools view -c "known_sites_chr4_chr16.vcf.gz" | bcftools reheader -h "contig_rename.txt" -O z -o "known_sites_chr4_chr16_NC.vcf.gz"
bcftools index "known_sites_chr4_chr16_NC.vcf.gz"

echo "Sorting & indexing BAM..."
samtools sort "SRR21384731_sorted_dedup_reads.bam" -o "SRR21384731_sorted_dedup_reads_sorted.bam"
samtools index "SRR21384731_sorted_dedup_reads_sorted.bam"

echo "Base Quality Score Recalibration - Table..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk BaseRecalibrator \
  -I /data/SRR21384731_sorted_dedup_reads_sorted.bam \
  --known-sites /data/known_sites_chr4_chr16_NC.vcf.gz \
  -O /data/SRR21384731_recal_data.table \
  -R /data/chr4_chr16.fa

echo "Base Quality Score Recalibration - Apply..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk ApplyBQSR \
  -R /data/chr4_chr16.fa \
  -I /data/SRR21384731_sorted_dedup_reads_sorted.bam \
  --bqsr-recal-file /data/SRR21384731_recal_data.table \
  -O /data/SRR21384731_sorted_dedup_BQSR.bam

echo "Calling variants with HaplotypeCaller..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk HaplotypeCaller \
  -R /data/chr4_chr16.fa \
  -I /data/SRR21384731_sorted_dedup_BQSR.bam \
  -O /data/SRR21384731_raw_variants.vcf \
  --standard-min-confidence-threshold-for-calling 20

echo "Split into SNPs and INDELs..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk SelectVariants \
  -R /data/chr4_chr16.fa \
  -V /data/SRR21384731_raw_variants.vcf \
  --select-type-to-include SNP \
  -O /data/SRR21384731_raw_snps.vcf

sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk SelectVariants \
  -R /data/chr4_chr16.fa \
  -V /data/SRR21384731_raw_variants.vcf \
  --select-type-to-include INDEL \
  -O /data/SRR21384731_raw_indels.vcf

echo "Hard-filter SNPs..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk VariantFiltration \
  -V /data/SRR21384731_raw_snps.vcf \
  --filter-expression "QD < 2.0" --filter-name "QD_filter" \
  --filter-expression "FS > 60.0" --filter-name "FS_filter" \
  --filter-expression "MQ < 40.0" --filter-name "MQ_filter" \
  --filter-expression "SOR > 4.0" --filter-name "SOR_filter" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_filter" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_filter" \
  --genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
  --genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter" \
  -O /data/SRR21384731_filtered_snps.vcf

echo "Selecting variants that pass filters with GATK..."
sudo docker run --rm -v "$PWD":/data "broadinstitute/gatk:latest" \
  gatk SelectVariants \
  --exclude-filtered \
  -V /data/SRR21384731_filtered_snps.vcf \
  -O /data/SRR21384731_analysis_ready_snps.vcf

echo "Converting VCF to ANNOVAR format..."
sudo docker run -it --rm -v "$PWD":/data "bioinfochrustrasbourg/annovar:latest" \
  perl ./convert2annovar.pl \
  -format vcf4 /data/SRR21384731_filtered_snps.vcf \
  -outfile /data/SRR21384731_filtered_snps.avinput

echo "Downloading ANNOVAR reference files..."
mkdir -p "humandb"
cd "humandb"
wget http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz
wget http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz

echo "ADPKD WES GATK Pipeline Completed Successfully!"
