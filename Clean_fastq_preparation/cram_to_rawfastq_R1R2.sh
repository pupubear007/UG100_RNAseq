#!/bin/bash

#SBATCH --job-name=cram2fastq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH --tmp=100gb
#SBATCH -o logs/cram2fastq_%j.out
#SBATCH -e logs/cram2fastq_%j.err

# CRAM to paired FASTQ conversion script
# Converts 120 unaligned CRAM files to R1/R2 FASTQ pairs
# Splits reads based on _1- (R1) and _2- (R2) in read names

# Load required modules
module load samtools
module load parallel

# Set paths
BASE_DIR="/scratch.global/wan00965/Quick_test"
CRAM_DIR="${BASE_DIR}/data/seqs"
FASTQ_DIR="${BASE_DIR}/data/fastq"
LOG_DIR="${BASE_DIR}/logs"
ANALYSIS_DIR="${BASE_DIR}/analysis"
SAMPLE_LOG="${ANALYSIS_DIR}/sample_log.csv"

# Set up TMPDIR in scratch (can delete manually later)
export TMPDIR="${BASE_DIR}/temp"
mkdir -p "${TMPDIR}"

echo "=== Environment Info ==="
echo "Temporary directory: ${TMPDIR}"
echo "Available temp space:"
df -h "${TMPDIR}" | tail -1
echo ""

# Create output directories
echo "=== Setting up directories ==="
mkdir -p "${FASTQ_DIR}"
mkdir -p "${LOG_DIR}/cram2fastq_conversion"

echo "CRAM directory: ${CRAM_DIR}"
echo "FASTQ output: ${FASTQ_DIR}"
echo "Temp directory: ${TMPDIR}"
echo "Sample metadata: ${SAMPLE_LOG}"
echo ""

# Check if sample log exists
if [ ! -f "${SAMPLE_LOG}" ]; then
    echo "WARNING: Sample log not found at ${SAMPLE_LOG}"
fi

# Change to CRAM directory
cd "${CRAM_DIR}" || exit 1

# Count total CRAM files
TOTAL_FILES=$(ls -1 *.cram 2>/dev/null | wc -l)
if [ ${TOTAL_FILES} -eq 0 ]; then
    echo "ERROR: No CRAM files found in ${CRAM_DIR}"
    exit 1
fi

echo "Found ${TOTAL_FILES} CRAM files to process"
echo ""

# Function to convert a single CRAM file to paired FASTQ
convert_cram_to_fastq() {
    local cram_file=$1
    local fastq_dir=$2
    local log_dir=$3
    local tmp_dir=$4
    
    local base_name=$(basename "${cram_file}" .cram)
    local log_file="${log_dir}/cram2fastq_conversion/${base_name}_conversion.log"
    
    # Start conversion
    {
        echo "=== Processing: ${base_name} ==="
        echo "Start time: $(date)"
        echo "Temp directory: ${tmp_dir}"
        echo ""
        
        # Check if output files already exist
        if [ -f "${fastq_dir}/${base_name}_R1.fastq.gz" ] && \
           [ -f "${fastq_dir}/${base_name}_R2.fastq.gz" ]; then
            echo "SKIPPED: FASTQ files already exist for ${base_name}"
            echo "End time: $(date)"
            return 0
        fi
        
        # Step 1: Convert CRAM to compressed interleaved FASTQ (saves space!)
        echo "Step 1: Extracting reads from CRAM to compressed temp file..."
        local temp_fastq="${tmp_dir}/${base_name}_all.fastq.gz"
        
        samtools fastq -@ 4 "${CRAM_DIR}/${cram_file}" | pigz -p 2 > "${temp_fastq}"
        
        local exit_code1=$?
        if [ ${exit_code1} -ne 0 ]; then
            echo "ERROR: Failed to extract FASTQ from CRAM (exit code: ${exit_code1})"
            rm -f "${temp_fastq}"
            return 1
        fi
        
        # Check if temp file was created and has content
        if [ ! -s "${temp_fastq}" ]; then
            echo "ERROR: Temp FASTQ file is empty or doesn't exist"
            return 1
        fi
        
        # Get file info
        local temp_size=$(du -h "${temp_fastq}" | cut -f1)
        echo "Temp file created: ${temp_fastq}"
        echo "Temp file size: ${temp_size}"
        
        # Count reads in compressed file
        local total_lines=$(zcat "${temp_fastq}" | wc -l)
        local total_reads=$((total_lines / 4))
        echo "Total reads extracted: ${total_reads}"
        echo ""
        
        # Step 2: Split into R1 and R2 based on _1- and _2- in read names
        echo "Step 2: Splitting into R1 and R2..."
        
        # Extract R1 (reads with _1- in header)
        zcat "${temp_fastq}" | \
        awk 'BEGIN {OFS="\n"} 
             /^@/ {header=$0; getline seq; getline plus; getline qual; 
                   if (header ~ /_1-/) print header,seq,plus,qual}' | \
        pigz -p 2 > "${fastq_dir}/${base_name}_R1.fastq.gz" &
        
        # Extract R2 (reads with _2- in header)
        zcat "${temp_fastq}" | \
        awk 'BEGIN {OFS="\n"} 
             /^@/ {header=$0; getline seq; getline plus; getline qual; 
                   if (header ~ /_2-/) print header,seq,plus,qual}' | \
        pigz -p 2 > "${fastq_dir}/${base_name}_R2.fastq.gz" &
        
        wait
        
        # Step 3: Verify outputs
        echo ""
        echo "Step 3: Verifying outputs..."
        
        local r1_lines=$(zcat "${fastq_dir}/${base_name}_R1.fastq.gz" | wc -l)
        local r2_lines=$(zcat "${fastq_dir}/${base_name}_R2.fastq.gz" | wc -l)
        local r1_reads=$((r1_lines / 4))
        local r2_reads=$((r2_lines / 4))
        
        echo "R1 reads: ${r1_reads}"
        echo "R2 reads: ${r2_reads}"
        echo "Total: $((r1_reads + r2_reads))"
        
        # Check if R1 and R2 counts are roughly equal
        if [ ${r1_reads} -eq 0 ] || [ ${r2_reads} -eq 0 ]; then
            echo "ERROR: One or both FASTQ files are empty!"
            echo "Keeping temp file for debugging: ${temp_fastq}"
            return 1
        fi
        
        local diff=$((r1_reads - r2_reads))
        if [ ${diff#-} -gt 1000 ]; then
            echo "WARNING: R1 and R2 read counts differ by ${diff}"
        fi
        
        # Get file sizes
        local r1_size=$(du -h "${fastq_dir}/${base_name}_R1.fastq.gz" | cut -f1)
        local r2_size=$(du -h "${fastq_dir}/${base_name}_R2.fastq.gz" | cut -f1)
        
        echo ""
        echo "SUCCESS: Converted ${base_name}"
        echo "  R1: ${r1_reads} reads (${r1_size})"
        echo "  R2: ${r2_reads} reads (${r2_size})"
        
        # Clean up temp file only if successful
        echo "  Removing temp file: ${temp_fastq}"
        rm -f "${temp_fastq}"
        
        echo "End time: $(date)"
        echo ""
        
        return 0
        
    } > "${log_file}" 2>&1
    
    local final_exit=$?
    
    # Print summary to stdout
    if [ ${final_exit} -eq 0 ]; then
        echo "✓ ${base_name}"
    else
        echo "✗ ${base_name} - FAILED"
    fi
    
    return ${final_exit}
}

# Export function and variables for GNU parallel
export -f convert_cram_to_fastq
export FASTQ_DIR LOG_DIR TMPDIR CRAM_DIR

# Process files in parallel
echo "=== Starting parallel CRAM to FASTQ conversion ==="
echo "Processing: $(date)"
echo ""

# With 4 threads per conversion and 16 CPUs total, run 4 conversions in parallel
PARALLEL_JOBS=$((SLURM_CPUS_PER_TASK / 4))
echo "Running ${PARALLEL_JOBS} parallel conversions (4 threads each)"
echo ""

ls -1 *.cram | parallel -j ${PARALLEL_JOBS} --progress \
    convert_cram_to_fastq {} "${FASTQ_DIR}" "${LOG_DIR}" "${TMPDIR}"

PARALLEL_EXIT=$?

echo ""
echo "=== Conversion complete ==="
echo "Finished: $(date)"
echo ""

# Generate summary statistics
echo "=== Summary Statistics ==="
echo ""

TOTAL_R1=$(ls -1 "${FASTQ_DIR}"/*_R1.fastq.gz 2>/dev/null | wc -l)
TOTAL_R2=$(ls -1 "${FASTQ_DIR}"/*_R2.fastq.gz 2>/dev/null | wc -l)

echo "Total CRAM files: ${TOTAL_FILES}"
echo "R1 files created: ${TOTAL_R1}"
echo "R2 files created: ${TOTAL_R2}"
echo "Complete pairs: $([ ${TOTAL_R1} -eq ${TOTAL_R2} ] && echo ${TOTAL_R1} || echo 'MISMATCH!')"
echo ""

# Check total size
echo "=== File size summary ==="
echo "FASTQ output directory:"
du -sh "${FASTQ_DIR}"
echo ""
echo "Temp directory (delete manually when done):"
du -sh "${TMPDIR}"
echo ""

# Count by host
echo "=== Samples by host ==="
SOYBEAN_R1=$(ls -1 "${FASTQ_DIR}"/*_Gm_*_R1.fastq.gz 2>/dev/null | wc -l)
SUNFLOWER_R1=$(ls -1 "${FASTQ_DIR}"/*_Ha_*_R1.fastq.gz 2>/dev/null | wc -l)
echo "  Soybean (Gm): ${SOYBEAN_R1}"
echo "  Sunflower (Ha): ${SUNFLOWER_R1}"
echo ""

# Identify any failed conversions
echo "=== Checking for failures ==="
FAILED_COUNT=0

for cram_file in *.cram; do
    base_name=$(basename "${cram_file}" .cram)
    if [ ! -f "${FASTQ_DIR}/${base_name}_R1.fastq.gz" ] || \
       [ ! -f "${FASTQ_DIR}/${base_name}_R2.fastq.gz" ]; then
        echo "MISSING: ${base_name}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
done

if [ ${FAILED_COUNT} -eq 0 ]; then
    echo "✓ All files converted successfully!"
else
    echo "✗ ${FAILED_COUNT} files failed to convert"
    echo "Check individual log files in: ${LOG_DIR}/cram2fastq_conversion/"
fi

echo ""
echo "=== Conversion logs location ==="
echo "${LOG_DIR}/cram2fastq_conversion/"
echo ""
echo "=== Temporary files location ==="
echo "${TMPDIR}"
echo "NOTE: Temp files are kept for manual inspection/deletion"
echo "To delete: rm -rf ${TMPDIR}"
echo ""

if [ ${PARALLEL_EXIT} -eq 0 ] && [ ${FAILED_COUNT} -eq 0 ]; then
    echo "✅ CRAM to FASTQ conversion completed successfully!"
    echo ""
    echo "Next steps:"
    echo "  1. Quality control: FastQC/MultiQC on FASTQ files"
    echo "  2. Alignment: Map reads to combined references using STAR"
    echo "  3. Clean up temp files: rm -rf ${TMPDIR}"
    exit 0
else
    echo "⚠️ CRAM to FASTQ conversion completed with errors"
    echo "Temp files kept for debugging in: ${TMPDIR}"
    exit 1
fi