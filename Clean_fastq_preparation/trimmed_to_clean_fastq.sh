#!/bin/bash
#SBATCH --job-name=cleanup_only
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --array=0-7%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/cleanup_%A_%a.out
#SBATCH -e logs/cleanup_%A_%a.err

# ============================================================================
# READ CLEANUP-ONLY PIPELINE
# 
# Processes raw reads through quality control and adapter removal
# Does NOT perform alignment - use alignment_only script after this!
# 
# Steps:
#   1. Rename and Repair (strip -XXXXXX, add /1 /2, repair.sh)
#   2. Cutadapt (V2 basic adapters, quality trimming)
#   3. FastQC (quality control on cleaned reads)
#   4. MultiQC (aggregate QC report)
# ============================================================================

set -e

echo "=== Cleanup-Only Pipeline - Node ${SLURM_ARRAY_TASK_ID} ==="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo ""

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR="/scratch.global/wan00965/Quick_test"

# Input/Output directories
INPUT_DIR="${BASE_DIR}/data/trimmed_fastq"           # Original paired files
TEMP_DIR="${BASE_DIR}/data/pipeline_temp"            # Temporary renaming
REPAIRED_DIR="${BASE_DIR}/data/repaired_reads"       # After repair.sh
CLEANED_DIR="${BASE_DIR}/data/cleaned_reads"         # After cutadapt
QC_DIR="${BASE_DIR}/analysis/qc"                     # QC reports

# Logs
LOG_DIR="${BASE_DIR}/logs/cleanup"

# Adapter file
ADAPTER_FILE="${BASE_DIR}/adapters/v2_basic_adapters.fasta"

# Create directories (only on first node)
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    mkdir -p "${TEMP_DIR}" "${REPAIRED_DIR}" "${CLEANED_DIR}"
    mkdir -p "${QC_DIR}"/{fastqc,multiqc_cleanup}
    mkdir -p "${LOG_DIR}"/{rename,cutadapt}
    mkdir -p logs
fi

# Load modules
module load conda
source activate pu_RNAseq

echo "✓ Environment loaded"
echo "✓ Cutadapt: $(cutadapt --version)"
echo "✓ FastQC: $(fastqc --version)"
echo ""

# Create adapter file if needed
if [ ! -f "${ADAPTER_FILE}" ]; then
    mkdir -p "$(dirname ${ADAPTER_FILE})"
    cat > ${ADAPTER_FILE} << 'EOF'
>overrep_seq_1
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCAGATCGGAAGAGCGTCG
>overrep_seq_2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGATCGGAAGAGCGTCG
>overrep_seq_3
GTGACTGGAGTTCAGACGTGTGCTCCTTCCGATCAGATCGGAAGAGCGTC
>overrep_seq_4
GTGACTGGAGTTCAGACGTGGTGCTCTTCCGATCAGATCGGAAGAGCGTC
>overrep_seq_5
GTGACTGGAGTTCAGACGTGGTGCTCCTTCCGATCAGATCGGAAGAGCGT
>overrep_seq_6
GTGACTGGAGTTCAGACGTGGTGCTCTTCCGATCTGATCGGAAGAGCGTC
>overrep_seq_7
GTGACTGGAGTTCAGACGTGTGCTCCTTCCGATCTGATCGGAAGAGCGTC
>overrep_seq_8
GCAGCAAGCAACTCTATCTCTCTAGAAAGGGGAGTGAGGGCCGAAAAATG
>overrep_seq_9
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTTTGTTGTGAGCTAT
>overrep_seq_10
GTGACTGAGTTCAGACGTGTGCTCTTCCGATCAGATCGGAAGAGCGTCGT
>overrep_seq_11
GTGACTGGAGTTCAGACGTGGTGCTCCTTCCGATCTGATCGGAAGAGCGT
>overrep_seq_12
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCTGCCAGTAGTCATA
>overrep_seq_13
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGAGTTCAGACGTGTG
>overrep_seq_14
GTGACTGAGTTCAGACGTGTGCTCTTCCGATCTGATCGGAAGAGCGTCGT
>Illumina_universal_adapter
AGATCGGAAGAG
EOF
    echo "✓ Created V2 adapter file"
fi

# ============================================================================
# Get samples for this node
# ============================================================================
ALL_SAMPLES=($(ls -1 "${INPUT_DIR}"/*_R1_paired.fastq.gz | sort))
TOTAL_SAMPLES=${#ALL_SAMPLES[@]}

if [ ${TOTAL_SAMPLES} -eq 0 ]; then
    echo "✗ ERROR: No samples found in ${INPUT_DIR}"
    echo "  Looking for: *_R1_paired.fastq.gz"
    exit 1
fi

echo "Found ${TOTAL_SAMPLES} samples in ${INPUT_DIR}"

SAMPLES_PER_NODE=$(((TOTAL_SAMPLES + 7) / 8))
START_IDX=$((SLURM_ARRAY_TASK_ID * SAMPLES_PER_NODE))
END_IDX=$((START_IDX + SAMPLES_PER_NODE - 1))

if [ ${END_IDX} -ge ${TOTAL_SAMPLES} ]; then
    END_IDX=$((TOTAL_SAMPLES - 1))
fi

if [ ${START_IDX} -ge ${TOTAL_SAMPLES} ]; then
    echo "Node ${SLURM_ARRAY_TASK_ID} has no samples assigned"
    exit 0
fi

NODE_SAMPLES=("${ALL_SAMPLES[@]:$START_IDX:$((END_IDX - START_IDX + 1))}")

echo "Node ${SLURM_ARRAY_TASK_ID}: Processing samples ${START_IDX} to ${END_IDX} (${#NODE_SAMPLES[@]} samples)"
echo ""

# ============================================================================
# Function to clean one sample
# ============================================================================
cleanup_sample() {
    local r1_input=$1
    local sample_name=$(basename "${r1_input}" _R1_paired.fastq.gz)
    local r2_input="${INPUT_DIR}/${sample_name}_R2_paired.fastq.gz"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Cleaning: ${sample_name}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Define file paths
    local renamed_r1="${TEMP_DIR}/${sample_name}_R1_renamed.fastq.gz"
    local renamed_r2="${TEMP_DIR}/${sample_name}_R2_renamed.fastq.gz"
    local repaired_r1="${REPAIRED_DIR}/${sample_name}_R1_repaired.fastq.gz"
    local repaired_r2="${REPAIRED_DIR}/${sample_name}_R2_repaired.fastq.gz"
    local singletons="${REPAIRED_DIR}/${sample_name}_singletons.fastq.gz"
    local cleaned_r1="${CLEANED_DIR}/${sample_name}_R1_cleaned.fastq.gz"
    local cleaned_r2="${CLEANED_DIR}/${sample_name}_R2_cleaned.fastq.gz"
    
    # ========================================================================
    # STEP 1: RENAME AND REPAIR
    # ========================================================================
    if [ ! -f "${repaired_r1}" ] || [ ! -f "${repaired_r2}" ]; then
        echo "[1/4] Rename and Repair..."
        
        # Strip -XXXXXX and add /1
        zcat "${r1_input}" | \
            awk 'NR%4==1 {sub(/-[0-9]+$/,""); print $0"/1"} NR%4!=1' | \
            gzip > "${renamed_r1}"
        
        # Strip -XXXXXX and add /2
        zcat "${r2_input}" | \
            awk 'NR%4==1 {sub(/-[0-9]+$/,""); print $0"/2"} NR%4!=1' | \
            gzip > "${renamed_r2}"
        
        # Run repair.sh
        repair.sh \
            in1="${renamed_r1}" \
            in2="${renamed_r2}" \
            out1="${repaired_r1}" \
            out2="${repaired_r2}" \
            outs="${singletons}" \
            repair=t \
            2>&1 | tee "${LOG_DIR}/rename/${sample_name}_repair.log"
        
        # Clean up temp files
        rm -f "${renamed_r1}" "${renamed_r2}"
        
        echo "    ✓ Repair complete"
    else
        echo "[1/4] ✓ Already repaired"
    fi
    
    # ========================================================================
    # STEP 2: CUTADAPT (ADAPTER REMOVAL + QUALITY TRIMMING)
    # ========================================================================
    if [ ! -f "${cleaned_r1}" ] || [ ! -f "${cleaned_r2}" ]; then
        echo "[2/4] Cutadapt (V2 basic adapters)..."
        
        cutadapt \
            -a file:${ADAPTER_FILE} \
            -A file:${ADAPTER_FILE} \
            -b file:${ADAPTER_FILE} \
            -B file:${ADAPTER_FILE} \
            -o "${cleaned_r1}" \
            -p "${cleaned_r2}" \
            --minimum-length 36 \
            --quality-cutoff 20 \
            --times 2 \
            -j 4 \
            "${repaired_r1}" "${repaired_r2}" \
            > "${LOG_DIR}/cutadapt/${sample_name}_cutadapt.txt" 2>&1
        
        echo "    ✓ Cutadapt complete"
    else
        echo "[2/4] ✓ Already cleaned"
    fi
    
    # ========================================================================
    # STEP 3: FASTQC ON CLEANED READS
    # ========================================================================
    echo "[3/4] FastQC on cleaned reads..."
    
    fastqc \
        -o "${QC_DIR}/fastqc" \
        -t 2 \
        --quiet \
        "${cleaned_r1}" "${cleaned_r2}"
    
    echo "    ✓ FastQC complete"
    
    # ========================================================================
    # GET STATISTICS
    # ========================================================================
    local retention="NA"
    if [ -f "${LOG_DIR}/cutadapt/${sample_name}_cutadapt.txt" ]; then
        local reads_in=$(grep "Total read pairs processed:" "${LOG_DIR}/cutadapt/${sample_name}_cutadapt.txt" | awk '{print $NF}' | sed 's/,//g')
        local reads_out=$(grep "Pairs written" "${LOG_DIR}/cutadapt/${sample_name}_cutadapt.txt" | awk '{print $NF}' | sed 's/,//g' | head -1)
        if [ -n "${reads_in}" ] && [ -n "${reads_out}" ]; then
            retention=$(awk "BEGIN {printf \"%.2f\", (${reads_out}/${reads_in})*100}")
        fi
    fi
    
    # Count reads in cleaned files
    local cleaned_reads=$(zcat "${cleaned_r1}" | wc -l | awk '{print $1/4}')
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "✓ ${sample_name} COMPLETE"
    echo "  Retention:     ${retention}%"
    echo "  Cleaned reads: $(printf "%'d" ${cleaned_reads}) pairs"
    echo "  Output:        ${CLEANED_DIR}/"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    return 0
}

export -f cleanup_sample
export INPUT_DIR TEMP_DIR REPAIRED_DIR CLEANED_DIR QC_DIR LOG_DIR ADAPTER_FILE

# ============================================================================
# Process samples (2 samples in parallel, each using 8 threads)
# ============================================================================
echo "=== Processing ${#NODE_SAMPLES[@]} Samples ==="
echo ""

printf "%s\n" "${NODE_SAMPLES[@]}" | \
    parallel -j 2 --bar cleanup_sample {}

# ============================================================================
# Node Summary
# ============================================================================
echo ""
echo "=== Node ${SLURM_ARRAY_TASK_ID} Complete ==="

COMPLETED=0
for r1_file in "${NODE_SAMPLES[@]}"; do
    sample_name=$(basename "${r1_file}" _R1_paired.fastq.gz)
    cleaned="${CLEANED_DIR}/${sample_name}_R1_cleaned.fastq.gz"
    
    if [ -f "${cleaned}" ]; then
        ((COMPLETED++))
    fi
done

printf "%-20s %d\n" "Assigned:" ${#NODE_SAMPLES[@]}
printf "%-20s %d\n" "Completed:" ${COMPLETED}
echo ""




######### MultiQC fail and need to debug
# ============================================================================
# Run MultiQC (Step 4 of 4) - only on last node
# ============================================================================
if [ ${SLURM_ARRAY_TASK_ID} -eq 7 ]; then
    echo "=== Generating Final MultiQC Report ==="
    sleep 120  # Wait for other nodes
    
    echo "[4/4] Running MultiQC on all cleaned reads..."
    cd "${BASE_DIR}/analysis"
    
    multiqc \
        "${QC_DIR}/fastqc" \
        "${LOG_DIR}/cutadapt" \
        -o "${QC_DIR}/multiqc_cleanup" \
        -n "cleanup_multiqc_report" \
        --title "Read Cleanup QC - All 120 Samples" \
        --force -q 2>&1
    
    echo "  ✓ MultiQC report: ${QC_DIR}/multiqc_cleanup/cleanup_multiqc_report.html"
    
    # Generate summary
    SUMMARY="${BASE_DIR}/analysis/cleanup_summary.txt"
    
    echo "=== Read Cleanup Pipeline Summary ===" > "${SUMMARY}"
    echo "Generated: $(date)" >> "${SUMMARY}"
    echo "" >> "${SUMMARY}"
    echo "Pipeline Steps:" >> "${SUMMARY}"
    echo "  1. Rename and Repair (strip -XXXXXX, add /1 /2, repair.sh)" >> "${SUMMARY}"
    echo "  2. Cutadapt (V2 basic adapters, Q20 trimming)" >> "${SUMMARY}"
    echo "  3. FastQC (quality control)" >> "${SUMMARY}"
    echo "  4. MultiQC (aggregate report)" >> "${SUMMARY}"
    echo "" >> "${SUMMARY}"
    
    # Count completions
    TOTAL_CLEANED=$(ls -1 "${CLEANED_DIR}"/*_R1_cleaned.fastq.gz 2>/dev/null | wc -l)
    
    echo "Samples completed: ${TOTAL_CLEANED}/${TOTAL_SAMPLES}" >> "${SUMMARY}"
    echo "" >> "${SUMMARY}"
    
    # Calculate average retention
    if [ ${TOTAL_CLEANED} -gt 0 ]; then
        echo "Read retention statistics:" >> "${SUMMARY}"
        
        TOTAL_RETENTION=0
        COUNT=0
        for cutadapt_log in "${LOG_DIR}"/cutadapt/*_cutadapt.txt; do
            if [ -f "${cutadapt_log}" ]; then
                reads_in=$(grep "Total read pairs processed:" "${cutadapt_log}" | awk '{print $NF}' | sed 's/,//g')
                reads_out=$(grep "Pairs written" "${cutadapt_log}" | awk '{print $NF}' | sed 's/,//g' | head -1)
                
                if [ -n "${reads_in}" ] && [ -n "${reads_out}" ] && [ "${reads_in}" -gt 0 ]; then
                    retention=$(awk "BEGIN {printf \"%.2f\", (${reads_out}/${reads_in})*100}")
                    TOTAL_RETENTION=$(awk "BEGIN {print ${TOTAL_RETENTION} + ${retention}}")
                    ((COUNT++))
                fi
            fi
        done
        
        if [ ${COUNT} -gt 0 ]; then
            AVG_RETENTION=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_RETENTION}/${COUNT}}")
            echo "  Average retention: ${AVG_RETENTION}%" >> "${SUMMARY}"
        fi
    fi
    
    echo "" >> "${SUMMARY}"
    echo "Output locations:" >> "${SUMMARY}"
    echo "  Repaired reads:  ${REPAIRED_DIR}/" >> "${SUMMARY}"
    echo "  Cleaned reads:   ${CLEANED_DIR}/" >> "${SUMMARY}"
    echo "  FastQC reports:  ${QC_DIR}/fastqc/" >> "${SUMMARY}"
    echo "  MultiQC report:  ${QC_DIR}/multiqc_cleanup/cleanup_multiqc_report.html" >> "${SUMMARY}"
    echo "" >> "${SUMMARY}"
    
    echo "✓ Summary: ${SUMMARY}"
    echo ""
    echo "✅ CLEANUP PIPELINE FINISHED!"
    echo ""
    echo "📊 Outputs:"
    echo "  Repaired reads:  ${REPAIRED_DIR}/"
    echo "  Cleaned reads:   ${CLEANED_DIR}/ ← USE THESE FOR ALIGNMENT"
    echo ""
    echo "  QC Reports:"
    echo "    FastQC:         ${QC_DIR}/fastqc/"
    echo "    MultiQC:        ${QC_DIR}/multiqc_cleanup/cleanup_multiqc_report.html"
    echo ""
    echo "  Summary:         ${SUMMARY}"
    echo ""
    echo "📝 Next Step:"
    echo "  Run alignment pipeline on cleaned reads:"
    echo "    CLEANED_DIR=\"${CLEANED_DIR}\"  # Set this in alignment script"
    echo "    sbatch alignment_only_pipeline_combined_ref.sh"
fi

echo ""
echo "End time: $(date)"
ELAPSED=$SECONDS
echo "Runtime: $(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"

exit 0
