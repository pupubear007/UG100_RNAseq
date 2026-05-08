#!/bin/bash
#SBATCH --job-name=star_only
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --array=0-7%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/alignment_%A_%a.out
#SBATCH -e logs/alignment_%A_%a.err

set -e

echo "=== STAR-Only Alignment Pipeline - Node ${SLURM_ARRAY_TASK_ID} ==="
echo "Start time: $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo ""

BASE_DIR="/scratch.global/wan00965/Quick_test"

# Input directory
CLEANED_DIR="${BASE_DIR}/data/U_final_cleaned_reads_re"

# Output directories
STAR_OUTPUT="${BASE_DIR}/analysis/star_alignment_combined_ref"

# Logs and statistics
LOG_DIR="${BASE_DIR}/logs/alignment"
STATS_DIR="${BASE_DIR}/analysis/alignment_stats"

# Reference indices - COMBINED VERSION
STAR_BIN="${BASE_DIR}/software/bin/STAR"
SOYBEAN_STAR="${BASE_DIR}/jgi_indices/star/combined/soybean_ss"
SUNFLOWER_STAR="${BASE_DIR}/jgi_indices/star/combined/sunflower_ss"

# GTF annotation files
if [ -f "${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gtf" ]; then
    SOYBEAN_GTF="${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gtf"
    SOYBEAN_FORMAT="GTF"
else
    SOYBEAN_GTF="${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gff3"
    SOYBEAN_FORMAT="GFF"
fi

if [ -f "${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gtf" ]; then
    SUNFLOWER_GTF="${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gtf"
    SUNFLOWER_FORMAT="GTF"
else
    SUNFLOWER_GTF="${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gff3"
    SUNFLOWER_FORMAT="GFF"
fi

# Create directories (only on first node)
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    mkdir -p "${STAR_OUTPUT}"/{soybean_samples,sunflower_samples}
    mkdir -p "${LOG_DIR}"
    mkdir -p "${STATS_DIR}"
    mkdir -p logs
fi

# Load modules
module load conda
source activate pu_RNAseq
module load samtools

echo "✓ Environment loaded"
echo "✓ STAR: $(${STAR_BIN} --version)"
echo "✓ Samtools: $(samtools --version | head -1)"
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
# Function to align one sample with STAR
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
    echo "Aligning: ${sample_name}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Determine host type and select indices
    local host_type
    local star_index gtf_file annotation_format
    
    if [[ ${sample_name} == *"_Gm_"* ]]; then
        host_type="soybean"
        star_index="${SOYBEAN_STAR}"
        gtf_file="${SOYBEAN_GTF}"
        annotation_format="${SOYBEAN_FORMAT}"
    elif [[ ${sample_name} == *"_Ha_"* ]]; then
        host_type="sunflower"
        star_index="${SUNFLOWER_STAR}"
        gtf_file="${SUNFLOWER_GTF}"
        annotation_format="${SUNFLOWER_FORMAT}"
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
    # STAR ALIGNMENT
    # ========================================================================
    local star_out="${STAR_OUTPUT}/${host_type}_samples/${sample_name}_"
    
    echo "[STAR] Starting alignment..."
    
    ${STAR_BIN} \
        --runThreadN 16 \
        --genomeDir "${star_index}" \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 100000000000 \
        --readFilesIn "${r1_input}" "${r2_input}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${star_out}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts \
        --sjdbGTFfile "${gtf_file}" \
        --outFilterMatchNminOverLread 0.2 \
        --outFilterScoreMinOverLread 0.2 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        > "${LOG_DIR}/${sample_name}_star.log" 2>&1
    
    # Index BAM
    samtools index "${star_out}Aligned.sortedByCoord.out.bam"
    
    # Flagstat
    samtools flagstat "${star_out}Aligned.sortedByCoord.out.bam" \
        > "${STATS_DIR}/${sample_name}_star_flagstat.txt"
    
    echo "    ✓ STAR complete"
    
    # ========================================================================
    # COLLECT STATISTICS
    # ========================================================================
    local star_mapped="NA"
    
    if [ -f "${STATS_DIR}/${sample_name}_star_flagstat.txt" ]; then
        star_mapped=$(grep "mapped (" "${STATS_DIR}/${sample_name}_star_flagstat.txt" | head -1 | awk '{print $5}' | sed 's/[()%]//g')
    fi
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "✓ ${sample_name} (${host_type}) COMPLETE"
    echo "  STAR:          ${star_mapped}%"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    return 0
}

export -f align_sample
export CLEANED_DIR STAR_OUTPUT LOG_DIR STATS_DIR
export STAR_BIN SOYBEAN_STAR SUNFLOWER_STAR SOYBEAN_GTF SUNFLOWER_GTF
export SOYBEAN_FORMAT SUNFLOWER_FORMAT

# ============================================================================
# Process samples sequentially (1 at a time, using all 16 threads)
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
        bam="${STAR_OUTPUT}/soybean_samples/${sample_name}_Aligned.sortedByCoord.out.bam"
    else
        bam="${STAR_OUTPUT}/sunflower_samples/${sample_name}_Aligned.sortedByCoord.out.bam"
    fi
    
    if [ -f "${bam}" ]; then
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
