#!/bin/bash

#SBATCH --job-name=fix_fastq
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/fix_fastq_%j.out
#SBATCH -e logs/fix_fastq_%j.err

# Fix Corrupted FASTQ Files
# Repairs FASTQ format issues that cause STAR alignment failures

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║           FASTQ FILE REPAIR UTILITY                         ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Start time: $(date)"
echo ""

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR="/scratch.global/wan00965/Quick_test"
CLEANED_DIR="${BASE_DIR}/data/final_cleaned_reads_re"
REPAIRED_DIR="${BASE_DIR}/data/sp_fixed_clean_reads"
LOG_DIR="${BASE_DIR}/logs/fastq_repair"

# Problem samples
PROBLEM_SAMPLES=(
    "48_3_Ha_SSPotter"
)


# Auto-discover all samples from the cleaned reads directory
echo "Discovering samples in: ${CLEANED_DIR}"

# Find all R1 files and extract sample names
PROBLEM_SAMPLES=($(find "${CLEANED_DIR}" -name "*_R1.fastq.gz" -exec basename {} \; | sed 's/_R1.fastq.gz//'))

echo "Found ${#PROBLEM_SAMPLES[@]} samples:"
for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo "  - ${sample}"
done
echo ""

# Verify that corresponding R2 files exist
echo "Verifying paired files..."
for sample in "${PROBLEM_SAMPLES[@]}"; do
    R1_FILE="${CLEANED_DIR}/${sample}_R1.fastq.gz"
    R2_FILE="${CLEANED_DIR}/${sample}_R2.fastq.gz"
    
    if [[ ! -f "${R1_FILE}" ]] || [[ ! -f "${R2_FILE}" ]]; then
        echo "⚠️  Warning: Missing pair for ${sample}"
        echo "   R1: $(test -f "${R1_FILE}" && echo "✓" || echo "✗")"
        echo "   R2: $(test -f "${R2_FILE}" && echo "✓" || echo "✗")"
    fi
done
echo ""



mkdir -p "${REPAIRED_DIR}"
mkdir -p "${LOG_DIR}"

# Load environment
module load conda
source activate pu_RNAseq

echo "✓ Environment loaded"
echo ""

# ============================================================================
# STEP 1: VALIDATE FASTQ FILES
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "STEP 1: Validate FASTQ Files"
echo "═══════════════════════════════════════════════════════════════"
echo ""

validate_fastq() {
    local file=$1
    local base_name=$(basename "${file}")
    
    echo "Validating: ${base_name}"
    
    # Use reformat.sh to validate
    reformat.sh \
        in="${file}" \
        out=/dev/null \
        2>&1 | tee "${LOG_DIR}/${base_name}_validation.log"
    
    if [ $? -eq 0 ]; then
        echo "✓ ${base_name} - Valid"
        return 0
    else
        echo "✗ ${base_name} - INVALID"
        return 1
    fi
}

echo "Checking problem samples:"
for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo ""
    echo "Sample: ${sample}"
    validate_fastq "${CLEANED_DIR}/${sample}_R1.fastq.gz"
    validate_fastq "${CLEANED_DIR}/${sample}_R2.fastq.gz"
done

echo ""

# ============================================================================
# STEP 2: IDENTIFY PROBLEMATIC READS
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "STEP 2: Identify Problematic Reads"
echo "═══════════════════════════════════════════════════════════════"
echo ""

find_bad_reads() {
    local file=$1
    local base_name=$(basename "${file}" .fastq.gz)
    local report="${LOG_DIR}/${base_name}_bad_reads.txt"
    
    echo "Scanning: ${base_name}"
    
    # Use awk to find mismatched seq/qual lengths
    zcat "${file}" | awk '
        NR % 4 == 1 { header = $0 }
        NR % 4 == 2 { seq = $0; seq_len = length($0) }
        NR % 4 == 0 { 
            qual_len = length($0)
            if (seq_len != qual_len) {
                print header "\tSeq=" seq_len "\tQual=" qual_len
            }
        }
    ' > "${report}"
    
    BAD_COUNT=$(wc -l < "${report}")
    
    if [ ${BAD_COUNT} -gt 0 ]; then
        echo "  Found ${BAD_COUNT} problematic reads"
        echo "  Details: ${report}"
        head -5 "${report}"
    else
        echo "  ✓ No format errors detected"
    fi
}

for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo ""
    find_bad_reads "${CLEANED_DIR}/${sample}_R1.fastq.gz"
    find_bad_reads "${CLEANED_DIR}/${sample}_R2.fastq.gz"
done

echo ""

# ============================================================================
# STEP 3: REPAIR FASTQ FILES
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "STEP 3: Repair FASTQ Files"
echo "═══════════════════════════════════════════════════════════════"
echo ""

repair_fastq() {
    local sample=$1
    local r1_in="${CLEANED_DIR}/${sample}_R1.fastq.gz"
    local r2_in="${CLEANED_DIR}/${sample}_R2.fastq.gz"
    local r1_out="${REPAIRED_DIR}/${sample}_R1_repaired.fastq.gz"
    local r2_out="${REPAIRED_DIR}/${sample}_R2_repaired.fastq.gz"
    local log="${LOG_DIR}/${sample}_repair.log"
    
    echo "Repairing: ${sample}"
    
    # Method 1: Use reformat.sh with strict validation
    reformat.sh \
        in1="${r1_in}" \
        in2="${r2_in}" \
        out1="${r1_out}" \
        out2="${r2_out}" \
        fixjunk=t \
        tossjunk=t \
        verifypaired=t \
        verifyinterleaved=f \
        allowidenticalnames=t \
        maxns=10 \
        minlen=36 \
        maxlen=1000 \
        threads=${SLURM_CPUS_PER_TASK} \
        2>&1 | tee "${log}"
    
    if [ $? -eq 0 ]; then
        echo "✓ ${sample} repaired successfully"
        
        # Verify output
        R1_COUNT=$(($(zcat "${r1_out}" | wc -l) / 4))
        R2_COUNT=$(($(zcat "${r2_out}" | wc -l) / 4))
        
        echo "  Output: ${R1_COUNT} read pairs"
        
        if [ ${R1_COUNT} -ne ${R2_COUNT} ]; then
            echo "  ⚠️  R1/R2 counts don't match!"
        fi
        
        return 0
    else
        echo "✗ ${sample} repair failed"
        return 1
    fi
}

for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo ""
    repair_fastq "${sample}"
done

echo ""

# ============================================================================
# STEP 4: ALTERNATIVE - STRICT FILTERING
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "STEP 4: Alternative - Aggressive Filtering"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "If repair didn't work, trying aggressive filtering..."
echo ""

aggressive_filter() {
    local sample=$1
    local r1_in="${CLEANED_DIR}/${sample}_R1.fastq.gz"
    local r2_in="${CLEANED_DIR}/${sample}_R2.fastq.gz"
    local r1_out="${REPAIRED_DIR}/${sample}_R1_filtered.fastq.gz"
    local r2_out="${REPAIRED_DIR}/${sample}_R2_filtered.fastq.gz"
    local log="${LOG_DIR}/${sample}_filter.log"
    
    echo "Aggressive filtering: ${sample}"
    
    # Very strict filtering
    reformat.sh \
        in1="${r1_in}" \
        in2="${r2_in}" \
        out1="${r1_out}" \
        out2="${r2_out}" \
        fixjunk=t \
        tossjunk=t \
        verifypaired=t \
        maxns=5 \
        minlen=50 \
        maxlen=500 \
        qtrim=rl \
        trimq=20 \
        minavgquality=20 \
        chastityfilter=t \
        barcodefilter=f \
        threads=${SLURM_CPUS_PER_TASK} \
        2>&1 | tee "${log}"
    
    if [ $? -eq 0 ]; then
        R1_COUNT=$(($(zcat "${r1_out}" | wc -l) / 4))
        echo "✓ Filtered: ${R1_COUNT} read pairs retained"
        return 0
    else
        echo "✗ Filtering failed"
        return 1
    fi
}

for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo ""
    aggressive_filter "${sample}"
done

echo ""

# ============================================================================
# STEP 5: VALIDATE REPAIRED FILES
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "STEP 5: Validate Repaired Files"
echo "═══════════════════════════════════════════════════════════════"
echo ""

for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo "Sample: ${sample}"
    
    # Check both repaired versions
    for suffix in "repaired" "filtered"; do
        r1="${REPAIRED_DIR}/${sample}_R1_${suffix}.fastq.gz"
        
        if [ -f "${r1}" ]; then
            echo ""
            echo "  Validating ${suffix} version:"
            
            # Validate format
            reformat.sh \
                in="${r1}" \
                out=/dev/null \
                verifypaired=f \
                2>&1 | grep -E "Input:|Output:" | head -2
            
            # Count reads
            COUNT=$(($(zcat "${r1}" | wc -l) / 4))
            printf "  Read count: %'d\n" ${COUNT}
        fi
    done
    
    echo ""
done

# ============================================================================
# STEP 6: RECOMMENDATION
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "RECOMMENDATIONS"
echo "═══════════════════════════════════════════════════════════════"
echo ""

for sample in "${PROBLEM_SAMPLES[@]}"; do
    echo "Sample: ${sample}"
    
    REPAIRED="${REPAIRED_DIR}/${sample}_R1_repaired.fastq.gz"
    FILTERED="${REPAIRED_DIR}/${sample}_R1_filtered.fastq.gz"
    
    if [ -f "${REPAIRED}" ]; then
        REPAIRED_COUNT=$(($(zcat "${REPAIRED}" | wc -l) / 4))
    else
        REPAIRED_COUNT=0
    fi
    
    if [ -f "${FILTERED}" ]; then
        FILTERED_COUNT=$(($(zcat "${FILTERED}" | wc -l) / 4))
    else
        FILTERED_COUNT=0
    fi
    
    echo "  Repaired version: $(printf '%\047d' ${REPAIRED_COUNT}) reads"
    echo "  Filtered version: $(printf '%\047d' ${FILTERED_COUNT}) reads"
    
    if [ ${REPAIRED_COUNT} -gt 0 ]; then
        echo "  ✓ RECOMMEND: Use repaired version"
        echo "    ${REPAIRED}"
    elif [ ${FILTERED_COUNT} -gt 0 ]; then
        echo "  ✓ RECOMMEND: Use filtered version"
        echo "    ${FILTERED}"
    else
        echo "  ✗ Both repairs failed - investigate original files"
    fi
    
    echo ""
done

# ============================================================================
# STEP 7: UPDATE SYMLINKS (Optional)
# ============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "UPDATE FOR STAR ALIGNMENT"
echo "═══════════════════════════════════════════════════════════════"
echo ""

echo "To use repaired files in STAR alignment:"
echo ""

for sample in "${PROBLEM_SAMPLES[@]}"; do
    REPAIRED_R1="${REPAIRED_DIR}/${sample}_R1_repaired.fastq.gz"
    REPAIRED_R2="${REPAIRED_DIR}/${sample}_R2_repaired.fastq.gz"
    
    if [ -f "${REPAIRED_R1}" ]; then
        echo "# For ${sample}:"
        echo "ln -sf ${REPAIRED_R1} ${CLEANED_DIR}/${sample}_R1_FIXED.fastq.gz"
        echo "ln -sf ${REPAIRED_R2} ${CLEANED_DIR}/${sample}_R2_FIXED.fastq.gz"
        echo ""
    fi
done

echo "Then update your STAR script to use *_FIXED.fastq.gz files"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                    REPAIR COMPLETE                           ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

echo "Output directory: ${REPAIRED_DIR}"
echo ""
echo "Next steps:"
echo "  1. Review logs in: ${LOG_DIR}"
echo "  2. Validate repaired files worked"
echo "  3. Re-run STAR with repaired files"
echo "  4. If still failing, check original raw data"
echo ""

echo "End time: $(date)"
echo ""