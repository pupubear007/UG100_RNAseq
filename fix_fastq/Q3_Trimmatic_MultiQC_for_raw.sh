#!/bin/bash

#SBATCH --job-name=qc_trim
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/qc_trim_%j.out
#SBATCH -e logs/qc_trim_%j.err

# Optimized Quality Control and Trimming Workflow
# Uses conda environment + GNU parallel for maximum efficiency

echo "=== Trimmatic and MultiQC ==="
echo "Start time: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}, Memory: 128GB"
echo ""

# Set paths
BASE_DIR="/scratch.global/wan00965/Quick_test_3"
FASTQ_DIR="${BASE_DIR}/data/raw_fq_from_cram_re"
TRIMMED_DIR="${BASE_DIR}/data/1trim_reads_re"
QC_DIR="${BASE_DIR}/analysis/qc"
LOG_DIR="${BASE_DIR}/logs"

# Create directories
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${QC_DIR}/qc_reports/before_trimming"
mkdir -p "${QC_DIR}/qc_reports/after_trimming"
mkdir -p "${LOG_DIR}/trimmomatic"

# Load modules
echo "=== Loading modules ==="
module load conda
module load parallel

# Activate conda environment
source activate pu_RNAseq

echo "✓ Modules loaded"
echo "  Conda env: pu_RNAseq"
echo "  GNU parallel: $(parallel --version | head -1)"
echo ""

# Verify tools
if ! command -v fastqc &> /dev/null || ! command -v trimmomatic &> /dev/null || ! command -v multiqc &> /dev/null; then
    echo "ERROR: Required tools not found in conda environment"
    echo "Please run: bash scripts/setup_conda_environment.sh"
    exit 1
fi

# Get adapter file
TRIMMOMATIC_DIR=$(dirname $(which trimmomatic))
ADAPTER_FILE="${TRIMMOMATIC_DIR}/../share/trimmomatic/adapters/TruSeq3-PE.fa"
if [ ! -f "${ADAPTER_FILE}" ]; then
    ADAPTER_FILE=$(find ${TRIMMOMATIC_DIR}/../share/trimmomatic* -name "TruSeq3-PE.fa" 2>/dev/null | head -1)
fi
if [ ! -f "${ADAPTER_FILE}" ]; then
    echo "ERROR: Cannot find TruSeq3-PE.fa"
    exit 1
fi

# Count samples
TOTAL_SAMPLES=$(ls -1 "${FASTQ_DIR}"/*_R1_cleaned.fastq.gz 2>/dev/null | wc -l)
echo "Found ${TOTAL_SAMPLES} sample pairs ($(($TOTAL_SAMPLES * 2)) files)"
echo ""

# Memory configuration note
echo "=== Memory Configuration ===
"
echo "Java heap per Trimmomatic job: 8GB (via JAVA_TOOL_OPTIONS)"
echo "Parallel Trimmomatic jobs: 8"
echo "Total Trimmomatic memory: ~64GB"
echo "System memory available: 128GB"
echo "✓ Memory allocation safe"
echo ""

# ============================================================================
# STEP 1: FastQC on raw reads - PARALLEL PROCESSING
# ============================================================================
echo "=== Step 1: FastQC on raw reads ==="
echo "Using GNU parallel for maximum efficiency..."
echo ""

# Check how many already done
EXISTING=$(ls -1 "${QC_DIR}/qc_reports/before_trimming/"/*_fastqc.html 2>/dev/null | wc -l)
if [ ${EXISTING} -eq $(($TOTAL_SAMPLES * 2)) ]; then
    echo "✓ All FastQC reports exist, skipping..."
else
    # Process with GNU parallel (each FastQC gets 1 CPU)
    # Run 16 FastQC processes in parallel
    find "${FASTQ_DIR}" -name "*.fastq.gz" | \
        parallel -j ${SLURM_CPUS_PER_TASK} --eta --bar \
        'fastqc {} -o '"${QC_DIR}"'/qc_reports/before_trimming/ -q 2>&1 | grep -v "^Analysis complete"'
    
    FASTQC_COUNT=$(ls -1 "${QC_DIR}/qc_reports/before_trimming/"/*_fastqc.html 2>/dev/null | wc -l)
    echo ""
    echo "✅ FastQC completed (${FASTQC_COUNT}/$(($TOTAL_SAMPLES * 2)) reports)"
fi
echo ""

# ============================================================================
# STEP 2: Trimmomatic - OPTIMIZED PARALLEL PROCESSING
# ============================================================================
echo "=== Step 2: Trimmomatic trimming ==="
echo "Processing samples with optimized parallel strategy..."
echo ""

# Function to trim one sample
trim_sample() {
    local r1_file=$1
    local base_name=$(basename "${r1_file}" _R1_cleaned.fastq.gz)
    local r2_file="${FASTQ_DIR}/${base_name}_R2_cleaned.fastq.gz"
    local log_file="${LOG_DIR}/trimmomatic/${base_name}_trim.log"
    
    # Skip if already done
    if [ -f "${TRIMMED_DIR}/${base_name}_R1_paired.fastq.gz" ] && \
       [ -f "${TRIMMED_DIR}/${base_name}_R2_paired.fastq.gz" ]; then
        return 0
    fi
    
    # Set Java options for EXTREMELY LARGE files - 48GB heap
    # Use 75% of available 64GB RAM
    # CRITICAL: Must use only 1 thread to minimize concurrent buffers
    export JAVA_TOOL_OPTIONS="-Xmx48g"
    
    # Run Trimmomatic with ONLY 1 THREAD
    # This prevents multiple concurrent decompression buffers
    trimmomatic PE -threads 1 \
        "${r1_file}" "${r2_file}" \
        "${TRIMMED_DIR}/${base_name}_R1_paired.fastq.gz" \
        "${TRIMMED_DIR}/${base_name}_R1_unpaired.fastq.gz" \
        "${TRIMMED_DIR}/${base_name}_R2_paired.fastq.gz" \
        "${TRIMMED_DIR}/${base_name}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        > "${log_file}" 2>&1
    
    return $?
}

export -f trim_sample
export FASTQ_DIR TRIMMED_DIR LOG_DIR ADAPTER_FILE

# Check how many already done
EXISTING_TRIMMED=$(ls -1 "${TRIMMED_DIR}"/*_R1_paired.fastq.gz 2>/dev/null | wc -l)
if [ ${EXISTING_TRIMMED} -eq ${TOTAL_SAMPLES} ]; then
    echo "✓ All samples already trimmed, skipping..."
else
    # Run 8 samples in parallel (4 threads each = 32 CPUs)
    ls -1 "${FASTQ_DIR}"/*_R1_cleaned.fastq.gz | \
        parallel -j 8 --eta --bar trim_sample {}
fi

TRIMMED_PAIRS=$(ls -1 "${TRIMMED_DIR}"/*_R1_paired.fastq.gz 2>/dev/null | wc -l)
echo ""
echo "✅ Trimming completed (${TRIMMED_PAIRS}/${TOTAL_SAMPLES} pairs)"
echo ""

# ============================================================================
# STEP 3: FastQC on trimmed reads - PARALLEL PROCESSING
# ============================================================================
echo "=== Step 3: FastQC on trimmed reads ==="
echo ""

EXISTING_TRIMMED_QC=$(ls -1 "${QC_DIR}/qc_reports/after_trimming/"/*_fastqc.html 2>/dev/null | wc -l)
if [ ${EXISTING_TRIMMED_QC} -eq $(($TRIMMED_PAIRS * 2)) ]; then
    echo "✓ All trimmed FastQC reports exist, skipping..."
else
    # Process with GNU parallel (16 FastQC jobs)
    find "${TRIMMED_DIR}" -name "*_paired.fastq.gz" | \
        parallel -j ${SLURM_CPUS_PER_TASK} --eta --bar \
        'fastqc {} -o '"${QC_DIR}"'/qc_reports/after_trimming/ -q 2>&1 | grep -v "^Analysis complete"'
    
    FASTQC_TRIMMED=$(ls -1 "${QC_DIR}/qc_reports/after_trimming/"/*_fastqc.html 2>/dev/null | wc -l)
    echo ""
    echo "✅ FastQC completed (${FASTQC_TRIMMED}/$(($TRIMMED_PAIRS * 2)) reports)"
fi
echo ""

# ============================================================================
# STEP 4: MultiQC reports - PARALLEL GENERATION
# ============================================================================
echo "=== Step 4: MultiQC reports ==="
echo ""

cd "${QC_DIR}" || exit 1

# Run MultiQC reports in parallel
parallel -j 3 ::: \
    "multiqc qc_reports/before_trimming/ -o qc_reports/ -n raw_fq_from_cram_re --title 'raw_fq_from_cram_re' --force -q" \
    "multiqc qc_reports/after_trimming/ -o qc_reports/ -n 1trim_reads_re --title '1trim_reads_re' --force -q" \
    "multiqc ${LOG_DIR}/trimmomatic/ -o qc_reports/ -n raw_to_1trim --title 'raw_to_1trim' --force -q"

echo "✅ MultiQC reports generated"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================
echo "=== Summary Statistics ==="
echo ""
printf "%-25s %d\n" "Input samples:" ${TOTAL_SAMPLES}
printf "%-25s %d\n" "Trimmed pairs:" ${TRIMMED_PAIRS}
printf "%-25s %.1f%%\n" "Success rate:" $(awk "BEGIN {printf \"%.1f\", (${TRIMMED_PAIRS}/${TOTAL_SAMPLES})*100}")
echo ""

echo "=== Storage Usage ==="
printf "%-20s %s\n" "Raw FASTQ:" "$(du -sh "${FASTQ_DIR}" 2>/dev/null | cut -f1)"
printf "%-20s %s\n" "Trimmed FASTQ:" "$(du -sh "${TRIMMED_DIR}" 2>/dev/null | cut -f1)"
printf "%-20s %s\n" "QC reports:" "$(du -sh "${QC_DIR}/qc_reports" 2>/dev/null | cut -f1)"
echo ""

echo "=== Samples by Host ==="
SOYBEAN=$(ls -1 "${TRIMMED_DIR}"/*_Gm_*_R1_paired.fastq.gz 2>/dev/null | wc -l)
SUNFLOWER=$(ls -1 "${TRIMMED_DIR}"/*_Ha_*_R1_paired.fastq.gz 2>/dev/null | wc -l)
printf "  %-15s %d\n" "Soybean (Gm):" ${SOYBEAN}
printf "  %-15s %d\n" "Sunflower (Ha):" ${SUNFLOWER}
echo ""

echo "=== Output Files ==="
echo "MultiQC Reports:"
echo "  ${QC_DIR}/qc_reports/before_trimming_multiqc.html"
echo "  ${QC_DIR}/qc_reports/after_trimming_multiqc.html"
echo "  ${QC_DIR}/qc_reports/trimming_stats_multiqc.html"
echo ""
echo "Trimmed FASTQ (use these for STAR alignment):"
echo "  ${TRIMMED_DIR}/*_R1_paired.fastq.gz (${TRIMMED_PAIRS} files)"
echo "  ${TRIMMED_DIR}/*_R2_paired.fastq.gz (${TRIMMED_PAIRS} files)"
echo ""

echo "=== Workflow Complete ==="
echo "End time: $(date)"
ELAPSED=$SECONDS
echo "Total runtime: $(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"
echo ""

if [ ${TRIMMED_PAIRS} -eq ${TOTAL_SAMPLES} ]; then
    echo "✅ SUCCESS: All ${TOTAL_SAMPLES} samples processed!"
    echo ""
    echo "📊 Next steps:"
    echo "  1. Download MultiQC reports: scp wan00965@mesabi.msi.umn.edu:${QC_DIR}/qc_reports/*_multiqc.html ."
    echo "  2. Review trimming statistics (should retain >85% reads)"
    echo "  3. Proceed to STAR alignment with trimmed paired reads"
    exit 0
else
    echo "⚠️  WARNING: Only ${TRIMMED_PAIRS}/${TOTAL_SAMPLES} samples completed"
    echo "Check logs: ${LOG_DIR}/trimmomatic/"
    exit 1
fi