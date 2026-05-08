#!/bin/bash
#SBATCH --job-name=salmon_only
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --array=0-7%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/salmon_%A_%a.out
#SBATCH -e logs/salmon_%A_%a.err

# ============================================================================
# SALMON-ONLY QUANTIFICATION PIPELINE - COMBINED REFERENCE VERSION
# 
# Starts from CLEANED READS (skip repair/cutadapt steps)
# Uses COMBINED indices (Host + Sclerotinia sclerotiorum)
# 
# Steps:
#   1. Salmon quantification
#   2. Collect statistics
# ============================================================================

set -e

echo "=== Salmon-Only Quantification Pipeline - Node ${SLURM_ARRAY_TASK_ID} ==="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo ""

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR="/scratch.global/wan00965/Quick_test"

# Input directory - CHANGE THIS to your cleaned reads location
CLEANED_DIR="${BASE_DIR}/data/U_final_cleaned_reads_re"

# Output directories
SALMON_OUTPUT="${BASE_DIR}/analysis/salmon_quant_combined_ref"
QC_DIR="${BASE_DIR}/analysis/qc"

# Logs and statistics
LOG_DIR="${BASE_DIR}/logs/salmon_quant"
STATS_DIR="${BASE_DIR}/analysis/salmon_stats"

# Reference indices - COMBINED VERSION
SOYBEAN_SALMON="${BASE_DIR}/jgi_indices/salmon/combined/soybean_ss_index"
SUNFLOWER_SALMON="${BASE_DIR}/jgi_indices/salmon/combined/sunflower_ss_index"

# Transcript-to-gene mapping for Salmon
SOYBEAN_TX2GENE="${BASE_DIR}/jgi_processed/tx2gene/soybean_ss_tx2gene.tsv"
SUNFLOWER_TX2GENE="${BASE_DIR}/jgi_processed/tx2gene/sunflower_ss_tx2gene.tsv"

# Create directories (only on first node)
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    mkdir -p "${SALMON_OUTPUT}"/{soybean_samples,sunflower_samples}
    mkdir -p "${LOG_DIR}"
    mkdir -p "${STATS_DIR}"
    mkdir -p logs
fi

# Load modules
module load conda
source activate pu_RNAseq
module load salmon

echo "✓ Environment loaded"
echo "✓ Salmon: $(salmon --version)"
echo ""

# ============================================================================
# Get samples for this node
# ============================================================================
ALL_SAMPLES=($(ls -1 "${CLEANED_DIR}"/*_R1.fastq.gz | sort))
TOTAL_SAMPLES=${#ALL_SAMPLES[@]}

if [ ${TOTAL_SAMPLES} -eq 0 ]; then
    echo "✗ ERROR: No samples found in ${CLEANED_DIR}"
    echo "  Looking for: *_R1.fastq.gz"
    echo ""
    echo "Available files:"
    ls -1 "${CLEANED_DIR}" | head -10
    exit 1
fi

echo "Found ${TOTAL_SAMPLES} samples in ${CLEANED_DIR}"

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
# Function to quantify one sample with Salmon
# ============================================================================
align_sample() {
    local r1_input=$1
    
    # Extract sample name
    local sample_name=$(basename "${r1_input}" | sed 's/_R1.*fastq.gz//')
    local r2_input=$(echo "${r1_input}" | sed 's/_R1/_R2/')
    
    # Verify R2 exists
    if [ ! -f "${r2_input}" ]; then
        echo "  ✗ ERROR: R2 file not found: ${r2_input}"
        return 1
    fi
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Quantifying: ${sample_name}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Determine host type and select indices
    local host_type
    local salmon_index tx2gene_file
    
    if [[ ${sample_name} == *"_Gm_"* ]]; then
        host_type="soybean"
        salmon_index="${SOYBEAN_SALMON}"
        tx2gene_file="${SOYBEAN_TX2GENE}"
    elif [[ ${sample_name} == *"_Ha_"* ]]; then
        host_type="sunflower"
        salmon_index="${SUNFLOWER_SALMON}"
        tx2gene_file="${SUNFLOWER_TX2GENE}"
    else
        echo "  ✗ Cannot determine host type for ${sample_name}"
        echo "  Expected '_Gm_' or '_Ha_' in sample name"
        return 1
    fi
    
    echo "Host: ${host_type}"
    echo "R1: ${r1_input}"
    echo "R2: ${r2_input}"
    echo ""
    
    # ========================================================================
    # SALMON QUANTIFICATION
    # ========================================================================
    local salmon_out="${SALMON_OUTPUT}/${host_type}_samples/${sample_name}"
    
    echo "[Salmon] Starting quantification..."
    
    salmon quant \
        -i "${salmon_index}" \
        -l A \
        -1 "${r1_input}" \
        -2 "${r2_input}" \
        -o "${salmon_out}" \
        --validateMappings \
        --minScoreFraction 0.6 \
        --minAssignedFrags 1 \
        -p 16 \
        -g "${tx2gene_file}" \
        > "${LOG_DIR}/${sample_name}_salmon.log" 2>&1
    
    echo "    ✓ Salmon complete (transcript + gene counts)"
    
    # ========================================================================
    # COLLECT STATISTICS
    # ========================================================================
    local salmon_mapped="NA"
    
    if [ -f "${salmon_out}/logs/salmon_quant.log" ]; then
        salmon_mapped=$(grep "Mapping rate" "${salmon_out}/logs/salmon_quant.log" | awk '{print $NF}' | sed 's/%//')
        
        # Copy log to stats directory for easier access
        cp "${salmon_out}/logs/salmon_quant.log" "${STATS_DIR}/${sample_name}_salmon_quant.log"
    fi
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "✓ ${sample_name} (${host_type}) COMPLETE"
    echo "  Salmon:        ${salmon_mapped}%"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    return 0
}

export -f align_sample
export CLEANED_DIR SALMON_OUTPUT LOG_DIR STATS_DIR
export SOYBEAN_SALMON SUNFLOWER_SALMON SOYBEAN_TX2GENE SUNFLOWER_TX2GENE

# ============================================================================
# Process samples sequentially
# ============================================================================
echo "=== Processing ${#NODE_SAMPLES[@]} Samples ==="
echo ""

for r1_file in "${NODE_SAMPLES[@]}"; do
    align_sample "${r1_file}"
done

# ============================================================================
# Node Summary
# ============================================================================
echo ""
echo "=== Node ${SLURM_ARRAY_TASK_ID} Complete ==="

COMPLETED=0
for r1_file in "${NODE_SAMPLES[@]}"; do
    sample_name=$(basename "${r1_file}" | sed 's/_R1.*fastq.gz//')
    
    if [[ ${sample_name} == *"_Gm_"* ]]; then
        quant_dir="${SALMON_OUTPUT}/soybean_samples/${sample_name}"
    else
        quant_dir="${SALMON_OUTPUT}/sunflower_samples/${sample_name}"
    fi
    
    if [ -f "${quant_dir}/quant.sf" ]; then
        ((COMPLETED++))
    fi
done

printf "%-20s %d\n" "Assigned:" ${#NODE_SAMPLES[@]}
printf "%-20s %d\n" "Completed:" ${COMPLETED}
echo ""

echo "End time: $(date)"
ELAPSED=$SECONDS
echo "Runtime: $(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"

exit 0
