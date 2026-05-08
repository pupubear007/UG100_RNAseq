#!/bin/bash
#SBATCH --job-name=build_all_jgi_indices
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=agsmall
#SBATCH -o logs/build_all_indices_%j.out
#SBATCH -e logs/build_all_indices_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu

# ========================================================================
# BUILD ALL JGI INDICES - Individual + Combined for Infection Studies
# ========================================================================

set -e

# Activate conda environment - EXACT USER SETTINGS
source ~/.bashrc
module load conda
source activate pu_RNAseq
# USE CUSTOM STAR COMPILED WITH DEF_readSeqLengthMax=700 (instead of default 650)
# This handles UG100 reads up to 653bp without trimming
export PATH="/scratch.global/wan00965/Quick_test_2/software/bin:${PATH}"
module load samtools
module load parallel
module load hisat2
module load salmon

echo "✓ STAR version: $(STAR --version) (custom build, max read length: 700bp)"
echo "✓ STAR location: $(which STAR)"
echo "✓ Samtools version: $(samtools --version | head -1)"
echo "✓ GNU parallel version: $(parallel --version | head -1)"
echo ""

# Verify custom STAR exists
if ! command -v STAR &> /dev/null; then
    echo "✗ ERROR: Custom STAR not found at /scratch.global/wan00965/Quick_test_2/software/bin/STAR"
    echo "Please compile STAR with DEF_readSeqLengthMax=700 first."
    exit 1
fi

# Directories
BASE_DIR="/scratch.global/wan00965/Quick_test_2"
JGI_DOWNLOAD="${BASE_DIR}/JGI_ref_genome_n_transcripts"
WORK_DIR="${BASE_DIR}/jgi_processed"
INDEX_DIR="${BASE_DIR}/jgi_indices"

# Create directory structure
mkdir -p ${WORK_DIR}/{soybean,sunflower,sclerotinia,combined/{soybean_ss,sunflower_ss}}
mkdir -p ${INDEX_DIR}/{star,hisat2,salmon}/{individual,combined}
mkdir -p ${BASE_DIR}/logs

echo "=========================================="
echo "JGI Genome Processing & Index Building"
echo "=========================================="
echo "Building 15 indices total:"
echo "  • 9 individual (3 organisms × 3 tools)"
echo "  • 6 combined (2 combos × 3 tools)"
echo ""
echo "CPUs: ${SLURM_CPUS_PER_TASK}, Memory: 64GB"
echo "Started: $(date)"
echo "=========================================="

# ========================================================================
# STEP 1: EXTRACT AND ORGANIZE FILES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 1: Extracting and Organizing Files"
echo "=========================================="

# ────────────────────────────────────────────────────────────────────────
# SOYBEAN - Glycine max (Wm82.a6.v1 from PhytozomeV13)
# ────────────────────────────────────────────────────────────────────────
echo ""
echo "Processing SOYBEAN files..."
cd ${WORK_DIR}/soybean

GM_GENOME="${JGI_DOWNLOAD}/Gm_genome/Phytozome/PhytozomeV13/Gmax/Wm82.a6.v1/assembly/Gmax_880_v6.0.fa.gz"
GM_TRANSCRIPT="${JGI_DOWNLOAD}/Gm_genome/Phytozome/PhytozomeV13/Gmax/Wm82.a6.v1/annotation/Gmax_880_Wm82.a6.v1.transcript.fa.gz"
GM_GFF="${JGI_DOWNLOAD}/Gm_genome/Phytozome/PhytozomeV13/Gmax/Wm82.a6.v1/annotation/Gmax_880_Wm82.a6.v1.gene.gff3.gz"

echo "  Extracting genome..."
gunzip -c ${GM_GENOME} > Gmax_genome.fa
echo "  ✓ Genome: $(du -h Gmax_genome.fa | cut -f1)"

echo "  Extracting transcripts..."
gunzip -c ${GM_TRANSCRIPT} > Gmax_transcripts.fa
echo "  ✓ Transcripts: $(du -h Gmax_transcripts.fa | cut -f1)"

echo "  Extracting annotations..."
gunzip -c ${GM_GFF} > Gmax_genes.gff3
echo "  ✓ GFF3: $(du -h Gmax_genes.gff3 | cut -f1)"

# ────────────────────────────────────────────────────────────────────────
# SUNFLOWER - Helianthus annuus (r1.2 from PhytozomeV12)
# ────────────────────────────────────────────────────────────────────────
echo ""
echo "Processing SUNFLOWER files..."
cd ${WORK_DIR}/sunflower

HA_GENOME="${JGI_DOWNLOAD}/Ha_genome/Phytozome/PhytozomeV12/early_release/Hannuus_494_r1.2/assembly/Hannuus_494_r1.0.fa.gz"
HA_TRANSCRIPT="${JGI_DOWNLOAD}/Ha_genome/Phytozome/PhytozomeV12/early_release/Hannuus_494_r1.2/annotation/Hannuus_494_r1.2.transcript.fa.gz"
HA_GFF="${JGI_DOWNLOAD}/Ha_genome/Phytozome/PhytozomeV12/early_release/Hannuus_494_r1.2/annotation/Hannuus_494_r1.2.gene.gff3.gz"

echo "  Extracting genome..."
gunzip -c ${HA_GENOME} > Hannuus_genome.fa
echo "  ✓ Genome: $(du -h Hannuus_genome.fa | cut -f1)"

echo "  Extracting transcripts..."
gunzip -c ${HA_TRANSCRIPT} > Hannuus_transcripts.fa
echo "  ✓ Transcripts: $(du -h Hannuus_transcripts.fa | cut -f1)"

echo "  Extracting annotations..."
gunzip -c ${HA_GFF} > Hannuus_genes.gff3
echo "  ✓ GFF3: $(du -h Hannuus_genes.gff3 | cut -f1)"

# ────────────────────────────────────────────────────────────────────────
# S. SCLEROTIORUM (from MycoCosm)
# ────────────────────────────────────────────────────────────────────────
echo ""
echo "Processing S. SCLEROTIORUM files..."
cd ${WORK_DIR}/sclerotinia

SS_GENOME="${JGI_DOWNLOAD}/Ss_genome/Mycocosm/Assembly/Genome Assembly (unmasked)/Sclsc1_AssemblyScaffolds.fasta.gz"
SS_TRANSCRIPT="${JGI_DOWNLOAD}/Ss_genome/Mycocosm/Annotation/Filtered Models (\"best\")/Transcripts/Sclsc1_GeneCatalog_transcripts_20110903.nt.fasta.gz"
SS_GFF="${JGI_DOWNLOAD}/Ss_genome/Mycocosm/Annotation/Filtered Models (\"best\")/Genes/Sclsc1_GeneCatalog_genes_20110903.gff.gz"

echo "  Extracting genome..."
gunzip -c "${SS_GENOME}" > Sclsc_genome.fa
echo "  ✓ Genome: $(du -h Sclsc_genome.fa | cut -f1)"

echo "  Extracting transcripts..."
gunzip -c "${SS_TRANSCRIPT}" > Sclsc_transcripts.fa
echo "  ✓ Transcripts: $(du -h Sclsc_transcripts.fa | cut -f1)"

echo "  Extracting annotations..."
gunzip -c "${SS_GFF}" > Sclsc_genes.gff
echo "  ✓ GFF: $(du -h Sclsc_genes.gff | cut -f1)"

echo ""
echo "✓ All files extracted successfully"

# ========================================================================
# STEP 2: CREATE COMBINED REFERENCES FOR INFECTION STUDIES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 2: Creating Combined References"
echo "=========================================="

# ────────────────────────────────────────────────────────────────────────
# COMBINED: SOYBEAN + S. SCLEROTIORUM
# ────────────────────────────────────────────────────────────────────────
echo ""
echo "Creating SOYBEAN + S. SCLEROTIORUM combined reference..."
cd ${WORK_DIR}/combined/soybean_ss

echo "  Combining genomes with prefixes..."
# Add Gm_ prefix to soybean FASTA headers to match GFF chromosome names
sed 's/^>/>Gm_/' ${WORK_DIR}/soybean/Gmax_genome.fa > Gm_prefixed_genome.fa
# Add Ss_ prefix to Sclerotinia FASTA headers to match GFF chromosome names
sed 's/^>/>Ss_/' ${WORK_DIR}/sclerotinia/Sclsc_genome.fa > Ss_prefixed_genome.fa
# Combine prefixed genomes
cat Gm_prefixed_genome.fa Ss_prefixed_genome.fa > Gm_Ss_genome.fa
echo "  ✓ Combined genome: $(du -h Gm_Ss_genome.fa | cut -f1)"

echo "  Combining transcripts with prefixes..."
sed 's/^>/>Gm_/' ${WORK_DIR}/soybean/Gmax_transcripts.fa > Gm_prefixed_transcripts.fa
sed 's/^>/>Ss_/' ${WORK_DIR}/sclerotinia/Sclsc_transcripts.fa > Ss_prefixed_transcripts.fa
cat Gm_prefixed_transcripts.fa Ss_prefixed_transcripts.fa > Gm_Ss_transcripts.fa
echo "  ✓ Combined transcripts: $(du -h Gm_Ss_transcripts.fa | cut -f1)"

echo "  Combining annotations with prefixes..."
awk '{if ($0 !~ /^#/) {print "Gm_" $0} else {print $0}}' ${WORK_DIR}/soybean/Gmax_genes.gff3 > Gm_prefixed_genes.gff3
awk '{if ($0 !~ /^#/) {print "Ss_" $0} else {print $0}}' ${WORK_DIR}/sclerotinia/Sclsc_genes.gff > Ss_prefixed_genes.gff
cat Gm_prefixed_genes.gff3 Ss_prefixed_genes.gff > Gm_Ss_genes.gff3
echo "  ✓ Combined annotations: $(du -h Gm_Ss_genes.gff3 | cut -f1)"

# ────────────────────────────────────────────────────────────────────────
# COMBINED: SUNFLOWER + S. SCLEROTIORUM
# ────────────────────────────────────────────────────────────────────────
echo ""
echo "Creating SUNFLOWER + S. SCLEROTIORUM combined reference..."
cd ${WORK_DIR}/combined/sunflower_ss

echo "  Combining genomes with prefixes..."
# Add Ha_ prefix to sunflower FASTA headers to match GFF chromosome names
sed 's/^>/>Ha_/' ${WORK_DIR}/sunflower/Hannuus_genome.fa > Ha_prefixed_genome.fa
# Add Ss_ prefix to Sclerotinia FASTA headers to match GFF chromosome names
sed 's/^>/>Ss_/' ${WORK_DIR}/sclerotinia/Sclsc_genome.fa > Ss_prefixed_genome.fa
# Combine prefixed genomes
cat Ha_prefixed_genome.fa Ss_prefixed_genome.fa > Ha_Ss_genome.fa
echo "  ✓ Combined genome: $(du -h Ha_Ss_genome.fa | cut -f1)"

echo "  Combining transcripts with prefixes..."
sed 's/^>/>Ha_/' ${WORK_DIR}/sunflower/Hannuus_transcripts.fa > Ha_prefixed_transcripts.fa
sed 's/^>/>Ss_/' ${WORK_DIR}/sclerotinia/Sclsc_transcripts.fa > Ss_prefixed_transcripts.fa
cat Ha_prefixed_transcripts.fa Ss_prefixed_transcripts.fa > Ha_Ss_transcripts.fa
echo "  ✓ Combined transcripts: $(du -h Ha_Ss_transcripts.fa | cut -f1)"

echo "  Combining annotations with prefixes..."
awk '{if ($0 !~ /^#/) {print "Ha_" $0} else {print $0}}' ${WORK_DIR}/sunflower/Hannuus_genes.gff3 > Ha_prefixed_genes.gff3
awk '{if ($0 !~ /^#/) {print "Ss_" $0} else {print $0}}' ${WORK_DIR}/sclerotinia/Sclsc_genes.gff > Ss_prefixed_genes.gff
cat Ha_prefixed_genes.gff3 Ss_prefixed_genes.gff > Ha_Ss_genes.gff3
echo "  ✓ Combined annotations: $(du -h Ha_Ss_genes.gff3 | cut -f1)"

echo ""
echo "✓ All combined references created"

# ========================================================================
# STEP 3: BUILD INDIVIDUAL STAR INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 3: Building Individual STAR Indices"
echo "=========================================="

STAR_INDIV_DIR="${INDEX_DIR}/star/individual"
mkdir -p ${STAR_INDIV_DIR}/{soybean,sunflower,sclerotinia}

# Soybean STAR
echo ""
echo "Building SOYBEAN STAR index..."
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDIV_DIR}/soybean \
    --genomeFastaFiles ${WORK_DIR}/soybean/Gmax_genome.fa \
    --sjdbGTFfile ${WORK_DIR}/soybean/Gmax_genes.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 652
echo "  ✓ Soybean STAR index complete"

# Sunflower STAR
echo ""
echo "Building SUNFLOWER STAR index..."
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDIV_DIR}/sunflower \
    --genomeFastaFiles ${WORK_DIR}/sunflower/Hannuus_genome.fa \
    --sjdbGTFfile ${WORK_DIR}/sunflower/Hannuus_genes.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 652
echo "  ✓ Sunflower STAR index complete"

# S. sclerotiorum STAR
echo ""
echo "Building S. SCLEROTIORUM STAR index..."
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDIV_DIR}/sclerotinia \
    --genomeFastaFiles ${WORK_DIR}/sclerotinia/Sclsc_genome.fa \
    --sjdbGTFfile ${WORK_DIR}/sclerotinia/Sclsc_genes.gff \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --genomeSAindexNbases 10 \
    --sjdbOverhang 652
echo "  ✓ S. sclerotiorum STAR index complete"

# ========================================================================
# STEP 4: BUILD COMBINED STAR INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 4: Building Combined STAR Indices"
echo "=========================================="

STAR_COMB_DIR="${INDEX_DIR}/star/combined"
mkdir -p ${STAR_COMB_DIR}/{soybean_ss,sunflower_ss}

# Soybean + S. sclerotiorum STAR
echo ""
echo "Building SOYBEAN + S. SCLEROTIORUM STAR index..."
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${STAR_COMB_DIR}/soybean_ss \
    --genomeFastaFiles ${WORK_DIR}/combined/soybean_ss/Gm_Ss_genome.fa \
    --sjdbGTFfile ${WORK_DIR}/combined/soybean_ss/Gm_Ss_genes.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 652
echo "  ✓ Soybean + Ss STAR index complete"

# Sunflower + S. sclerotiorum STAR
echo ""
echo "Building SUNFLOWER + S. SCLEROTIORUM STAR index..."
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ${STAR_COMB_DIR}/sunflower_ss \
    --genomeFastaFiles ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_genome.fa \
    --sjdbGTFfile ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_genes.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon CDS \
    --genomeSAindexNbases 13 \
    --sjdbOverhang 652
echo "  ✓ Sunflower + Ss STAR index complete"

# ========================================================================
# STEP 5: BUILD INDIVIDUAL HISAT2 INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 5: Building Individual HISAT2 Indices"
echo "=========================================="

HISAT2_INDIV_DIR="${INDEX_DIR}/hisat2/individual"
mkdir -p ${HISAT2_INDIV_DIR}

# Soybean HISAT2
echo ""
echo "Building SOYBEAN HISAT2 index..."
cd ${HISAT2_INDIV_DIR}
hisat2_extract_splice_sites.py ${WORK_DIR}/soybean/Gmax_genes.gff3 > soybean_splicesites.txt
hisat2_extract_exons.py ${WORK_DIR}/soybean/Gmax_genes.gff3 > soybean_exons.txt
hisat2-build \
    -p 16 \
    --ss soybean_splicesites.txt \
    --exon soybean_exons.txt \
    ${WORK_DIR}/soybean/Gmax_genome.fa \
    soybean_index
echo "  ✓ Soybean HISAT2 index complete"

# Sunflower HISAT2
echo ""
echo "Building SUNFLOWER HISAT2 index..."
hisat2_extract_splice_sites.py ${WORK_DIR}/sunflower/Hannuus_genes.gff3 > sunflower_splicesites.txt
hisat2_extract_exons.py ${WORK_DIR}/sunflower/Hannuus_genes.gff3 > sunflower_exons.txt
hisat2-build \
    -p 16 \
    --ss sunflower_splicesites.txt \
    --exon sunflower_exons.txt \
    ${WORK_DIR}/sunflower/Hannuus_genome.fa \
    sunflower_index
echo "  ✓ Sunflower HISAT2 index complete"

# S. sclerotiorum HISAT2
echo ""
echo "Building S. SCLEROTIORUM HISAT2 index..."
hisat2_extract_splice_sites.py ${WORK_DIR}/sclerotinia/Sclsc_genes.gff > sclerotinia_splicesites.txt
hisat2_extract_exons.py ${WORK_DIR}/sclerotinia/Sclsc_genes.gff > sclerotinia_exons.txt
hisat2-build \
    -p 16 \
    --ss sclerotinia_splicesites.txt \
    --exon sclerotinia_exons.txt \
    ${WORK_DIR}/sclerotinia/Sclsc_genome.fa \
    sclerotinia_index
echo "  ✓ S. sclerotiorum HISAT2 index complete"

# ========================================================================
# STEP 6: BUILD COMBINED HISAT2 INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 6: Building Combined HISAT2 Indices"
echo "=========================================="

HISAT2_COMB_DIR="${INDEX_DIR}/hisat2/combined"
mkdir -p ${HISAT2_COMB_DIR}

# Soybean + S. sclerotiorum HISAT2
echo ""
echo "Building SOYBEAN + S. SCLEROTIORUM HISAT2 index..."
cd ${HISAT2_COMB_DIR}
hisat2_extract_splice_sites.py ${WORK_DIR}/combined/soybean_ss/Gm_Ss_genes.gff3 > soybean_ss_splicesites.txt
hisat2_extract_exons.py ${WORK_DIR}/combined/soybean_ss/Gm_Ss_genes.gff3 > soybean_ss_exons.txt
hisat2-build \
    -p 16 \
    --ss soybean_ss_splicesites.txt \
    --exon soybean_ss_exons.txt \
    ${WORK_DIR}/combined/soybean_ss/Gm_Ss_genome.fa \
    soybean_ss_index
echo "  ✓ Soybean + Ss HISAT2 index complete"

# Sunflower + S. sclerotiorum HISAT2
echo ""
echo "Building SUNFLOWER + S. SCLEROTIORUM HISAT2 index..."
hisat2_extract_splice_sites.py ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_genes.gff3 > sunflower_ss_splicesites.txt
hisat2_extract_exons.py ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_genes.gff3 > sunflower_ss_exons.txt
hisat2-build \
    -p 16 \
    --ss sunflower_ss_splicesites.txt \
    --exon sunflower_ss_exons.txt \
    ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_genome.fa \
    sunflower_ss_index
echo "  ✓ Sunflower + Ss HISAT2 index complete"

# ========================================================================
# STEP 7: BUILD INDIVIDUAL SALMON INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 7: Building Individual Salmon Indices"
echo "=========================================="

SALMON_INDIV_DIR="${INDEX_DIR}/salmon/individual"
mkdir -p ${SALMON_INDIV_DIR}

# Soybean Salmon
echo ""
echo "Building SOYBEAN Salmon index..."
salmon index \
    -t ${WORK_DIR}/soybean/Gmax_transcripts.fa \
    -i ${SALMON_INDIV_DIR}/soybean_index \
    -p 16
echo "  ✓ Soybean Salmon index complete"

# Sunflower Salmon
echo ""
echo "Building SUNFLOWER Salmon index..."
salmon index \
    -t ${WORK_DIR}/sunflower/Hannuus_transcripts.fa \
    -i ${SALMON_INDIV_DIR}/sunflower_index \
    -p 16
echo "  ✓ Sunflower Salmon index complete"

# S. sclerotiorum Salmon
echo ""
echo "Building S. SCLEROTIORUM Salmon index..."
salmon index \
    -t ${WORK_DIR}/sclerotinia/Sclsc_transcripts.fa \
    -i ${SALMON_INDIV_DIR}/sclerotinia_index \
    -p 16
echo "  ✓ S. sclerotiorum Salmon index complete"

# ========================================================================
# STEP 8: BUILD COMBINED SALMON INDICES
# ========================================================================
echo ""
echo "=========================================="
echo "STEP 8: Building Combined Salmon Indices"
echo "=========================================="

SALMON_COMB_DIR="${INDEX_DIR}/salmon/combined"
mkdir -p ${SALMON_COMB_DIR}

# Soybean + S. sclerotiorum Salmon
echo ""
echo "Building SOYBEAN + S. SCLEROTIORUM Salmon index..."
salmon index \
    -t ${WORK_DIR}/combined/soybean_ss/Gm_Ss_transcripts.fa \
    -i ${SALMON_COMB_DIR}/soybean_ss_index \
    -p 16
echo "  ✓ Soybean + Ss Salmon index complete"

# Sunflower + S. sclerotiorum Salmon
echo ""
echo "Building SUNFLOWER + S. SCLEROTIORUM Salmon index..."
salmon index \
    -t ${WORK_DIR}/combined/sunflower_ss/Ha_Ss_transcripts.fa \
    -i ${SALMON_COMB_DIR}/sunflower_ss_index \
    -p 16
echo "  ✓ Sunflower + Ss Salmon index complete"

# ========================================================================
# FINAL SUMMARY
# ========================================================================
echo ""
echo "=========================================="
echo "🎉 ALL INDICES BUILT SUCCESSFULLY!"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "──────────────────────────────────────────"
echo "INDIVIDUAL INDICES (9 total):"
echo "──────────────────────────────────────────"
echo ""
echo "STAR Individual:"
echo "  • ${STAR_INDIV_DIR}/soybean/"
echo "  • ${STAR_INDIV_DIR}/sunflower/"
echo "  • ${STAR_INDIV_DIR}/sclerotinia/"
echo ""
echo "HISAT2 Individual:"
echo "  • ${HISAT2_INDIV_DIR}/soybean_index"
echo "  • ${HISAT2_INDIV_DIR}/sunflower_index"
echo "  • ${HISAT2_INDIV_DIR}/sclerotinia_index"
echo ""
echo "Salmon Individual:"
echo "  • ${SALMON_INDIV_DIR}/soybean_index/"
echo "  • ${SALMON_INDIV_DIR}/sunflower_index/"
echo "  • ${SALMON_INDIV_DIR}/sclerotinia_index/"
echo ""
echo "──────────────────────────────────────────"
echo "COMBINED INDICES (6 total):"
echo "──────────────────────────────────────────"
echo ""
echo "STAR Combined:"
echo "  • ${STAR_COMB_DIR}/soybean_ss/"
echo "  • ${STAR_COMB_DIR}/sunflower_ss/"
echo ""
echo "HISAT2 Combined:"
echo "  • ${HISAT2_COMB_DIR}/soybean_ss_index"
echo "  • ${HISAT2_COMB_DIR}/sunflower_ss_index"
echo ""
echo "Salmon Combined:"
echo "  • ${SALMON_COMB_DIR}/soybean_ss_index/"
echo "  • ${SALMON_COMB_DIR}/sunflower_ss_index/"
echo ""
echo "──────────────────────────────────────────"
echo "USAGE GUIDE:"
echo "──────────────────────────────────────────"
echo ""
echo "For CONTROL samples (0 hpi - no infection):"
echo "  → Use INDIVIDUAL indices"
echo ""
echo "For INFECTED samples (24, 48, 96 hpi):"
echo "  → Use COMBINED indices"
echo ""
echo "──────────────────────────────────────────"
echo "NEXT STEP: Test alignment"
echo "  sbatch test_all_jgi_indices.sh"
echo "=========================================="
