#!/bin/bash
#SBATCH --job-name=dual_count_v5
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --array=0-7%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/dual_count_v5_%A_%a.out
#SBATCH -e logs/dual_count_v5_%A_%a.err

set -e

echo "=== Dual Counting Pipeline v5 (correct attribute keys per file) - Node ${SLURM_ARRAY_TASK_ID} ==="
echo "Start time: $(date)"
echo ""

BASE_DIR="/scratch.global/wan00965/Quick_test"
STAR_OUTPUT="${BASE_DIR}/analysis/star_alignment_combined_ref"
SOYBEAN_BAM_DIR="${STAR_OUTPUT}/soybean_samples"
SUNFLOWER_BAM_DIR="${STAR_OUTPUT}/sunflower_samples"
COUNTS_DIR="${BASE_DIR}/analysis/dual_counts_v5"
LOG_DIR="${BASE_DIR}/logs/dual_counts_v5"
ANNO_DIR="${BASE_DIR}/jgi_processed/combined"

# ============================================================================
# Setup (only on first node)
# ============================================================================
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    mkdir -p "${COUNTS_DIR}"/{soybean,sunflower}
    mkdir -p "${ANNO_DIR}/split_annotations"
    mkdir -p "${LOG_DIR}"
    mkdir -p logs

    echo "Creating split annotation files..."

    if [ ! -f "${ANNO_DIR}/split_annotations/soybean_only.gff3" ]; then
        grep -E "^#|^Gm_" "${ANNO_DIR}/soybean_ss/Gm_Ss_genes.gff3" \
            > "${ANNO_DIR}/split_annotations/soybean_only.gff3"
        echo "  ✓ soybean_only.gff3"
    fi

    if [ ! -f "${ANNO_DIR}/split_annotations/pathogen_soy.gff3" ]; then
        grep -E "^#|^Ss_" "${ANNO_DIR}/soybean_ss/Gm_Ss_genes.gff3" \
            > "${ANNO_DIR}/split_annotations/pathogen_soy.gff3"
        echo "  ✓ pathogen_soy.gff3"
    fi

    if [ ! -f "${ANNO_DIR}/split_annotations/sunflower_only.gff3" ]; then
        grep -E "^#|^Ha_" "${ANNO_DIR}/sunflower_ss/Ha_Ss_genes.gff3" \
            > "${ANNO_DIR}/split_annotations/sunflower_only.gff3"
        echo "  ✓ sunflower_only.gff3"
    fi

    if [ ! -f "${ANNO_DIR}/split_annotations/pathogen_sun.gff3" ]; then
        grep -E "^#|^Ss_" "${ANNO_DIR}/sunflower_ss/Ha_Ss_genes.gff3" \
            > "${ANNO_DIR}/split_annotations/pathogen_sun.gff3"
        echo "  ✓ pathogen_sun.gff3"
    fi

    # Show one CDS line from each so we can verify attribute keys
    echo ""
    echo "=== Annotation attribute samples (verify before running) ==="
    for f in soybean_only.gff3 pathogen_soy.gff3 sunflower_only.gff3 pathogen_sun.gff3; do
        echo "--- ${f} ---"
        grep -P "\tCDS\t" "${ANNO_DIR}/split_annotations/${f}" 2>/dev/null \
            | head -1 | cut -f9 || echo "  (no CDS records found!)"
    done
    echo ""
fi

sleep 5

# ============================================================================
# Load environment
# ============================================================================
module load conda
source activate pu_RNAseq
module load parallel

echo "✓ Environment loaded"
echo "✓ featureCounts: $(featureCounts -v 2>&1 | head -1)"
echo ""

# ============================================================================
# Gather samples
# ============================================================================
ALL_SOYBEAN_BAMS=($(ls -1 "${SOYBEAN_BAM_DIR}"/*_Aligned.sortedByCoord.out.bam 2>/dev/null | sort))
ALL_SUNFLOWER_BAMS=($(ls -1 "${SUNFLOWER_BAM_DIR}"/*_Aligned.sortedByCoord.out.bam 2>/dev/null | sort))

TOTAL_SOYBEAN=${#ALL_SOYBEAN_BAMS[@]}
TOTAL_SUNFLOWER=${#ALL_SUNFLOWER_BAMS[@]}
TOTAL_SAMPLES=$((TOTAL_SOYBEAN + TOTAL_SUNFLOWER))

echo "Found ${TOTAL_SOYBEAN} soybean samples"
echo "Found ${TOTAL_SUNFLOWER} sunflower samples"
echo ""

ALL_SAMPLES=("${ALL_SOYBEAN_BAMS[@]}" "${ALL_SUNFLOWER_BAMS[@]}")

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

echo "Node ${SLURM_ARRAY_TASK_ID}: Processing ${#NODE_SAMPLES[@]} samples"
echo ""

# ============================================================================
# Function: count one sample
# Host  uses -g Parent  (GFF3-style attributes:  ID=...;Parent=...)
# Pathogen uses -g name (GTF-style attributes:   name "SS1G_..."; proteinId 1)
# Pathogen also uses -t exon (GTF from JGI typically lacks 'CDS' rows but has 'exon')
# ============================================================================
count_sample() {
    local bam_file=$1
    local sample_name=$(basename "${bam_file}" | sed 's/_Aligned.sortedByCoord.out.bam//')

    local host_anno pathogen_anno output_dir
    local host_feature host_attr pathogen_feature pathogen_attr

    if [[ ${bam_file} == *"/soybean_samples/"* ]]; then
        host_anno="${ANNO_DIR}/split_annotations/soybean_only.gff3"
        pathogen_anno="${ANNO_DIR}/split_annotations/pathogen_soy.gff3"
        output_dir="${COUNTS_DIR}/soybean"
    elif [[ ${bam_file} == *"/sunflower_samples/"* ]]; then
        host_anno="${ANNO_DIR}/split_annotations/sunflower_only.gff3"
        pathogen_anno="${ANNO_DIR}/split_annotations/pathogen_sun.gff3"
        output_dir="${COUNTS_DIR}/sunflower"
    else
        echo "[${sample_name}] ERROR: Cannot determine host type"
        return 1
    fi

    # Host annotation: GFF3 from JGI (Glyma/HanXRQ); CDS feature, Parent attribute
    host_feature="CDS"
    host_attr="Parent"

    # Pathogen annotation: GTF-style content from JGI Sclerotinia; uses 'name'
    # JGI fungal GTFs typically have 'exon' rows with 'name "SS1G_..."'
    pathogen_feature="exon"
    pathogen_attr="name"

    local final_output="${output_dir}/${sample_name}_combined_counts.txt"

    if [ -f "${final_output}" ]; then
        echo "[${sample_name}] Already processed"
        return 0
    fi

    echo "[${sample_name}] Starting..."
    local temp_dir="${output_dir}/temp_${sample_name}"
    mkdir -p "${temp_dir}"

    # ------------------------------------------------------------------------
    # Step 1: HOST
    # ------------------------------------------------------------------------
    echo "[${sample_name}] Counting host genes (-t ${host_feature} -g ${host_attr})..."
    featureCounts \
        -p \
        -B \
        -C \
        -T 4 \
        -t "${host_feature}" \
        -g "${host_attr}" \
        -a "${host_anno}" \
        -o "${temp_dir}/host_counts.txt" \
        "${bam_file}" \
        > "${LOG_DIR}/${sample_name}_host.log" 2>&1

    if [ $? -ne 0 ]; then
        echo "[${sample_name}] Host counting FAILED — see ${LOG_DIR}/${sample_name}_host.log"
        rm -rf "${temp_dir}"; return 1
    fi

    # ------------------------------------------------------------------------
    # Step 2: PATHOGEN
    # ------------------------------------------------------------------------
    echo "[${sample_name}] Counting pathogen genes (-t ${pathogen_feature} -g ${pathogen_attr})..."
    featureCounts \
        -p \
        -B \
        -C \
        -T 4 \
        -t "${pathogen_feature}" \
        -g "${pathogen_attr}" \
        -a "${pathogen_anno}" \
        -o "${temp_dir}/pathogen_counts.txt" \
        "${bam_file}" \
        > "${LOG_DIR}/${sample_name}_pathogen.log" 2>&1

    if [ $? -ne 0 ]; then
        echo "[${sample_name}] Pathogen counting FAILED — see ${LOG_DIR}/${sample_name}_pathogen.log"
        rm -rf "${temp_dir}"; return 1
    fi

    # ------------------------------------------------------------------------
    # Step 3: Merge
    # featureCounts output: line 1 = comment, line 2 = column header,
    # then GeneID Chr Start End Strand Length <bam_path> (count is col 7)
    # ------------------------------------------------------------------------
    awk 'NR>2 {print $1"\t"$7}' "${temp_dir}/host_counts.txt"     > "${temp_dir}/host_clean.txt"
    awk 'NR>2 {print $1"\t"$7}' "${temp_dir}/pathogen_counts.txt" > "${temp_dir}/pathogen_clean.txt"

    echo -e "GeneID\t${sample_name}" > "${final_output}"
    cat "${temp_dir}/host_clean.txt" "${temp_dir}/pathogen_clean.txt" >> "${final_output}"

    local HOST_GENES=$(wc -l < "${temp_dir}/host_clean.txt")
    local PATH_GENES=$(wc -l < "${temp_dir}/pathogen_clean.txt")
    local HOST_ASSIGNED=$(awk '$1=="Assigned"{print $2}' "${temp_dir}/host_counts.txt.summary" 2>/dev/null || echo "?")
    local PATH_ASSIGNED=$(awk '$1=="Assigned"{print $2}' "${temp_dir}/pathogen_counts.txt.summary" 2>/dev/null || echo "?")

    echo "[${sample_name}] ✓ Done: ${HOST_GENES} host + ${PATH_GENES} pathogen genes | reads: host=${HOST_ASSIGNED}, pathogen=${PATH_ASSIGNED}"

    cp "${temp_dir}/host_counts.txt.summary"     "${output_dir}/${sample_name}_host_summary.txt"     2>/dev/null || true
    cp "${temp_dir}/pathogen_counts.txt.summary" "${output_dir}/${sample_name}_pathogen_summary.txt" 2>/dev/null || true

    rm -rf "${temp_dir}"
    return 0
}

export -f count_sample
export COUNTS_DIR LOG_DIR ANNO_DIR BASE_DIR

# ============================================================================
# Run with parallel
# ============================================================================
echo "=== Starting Parallel Processing ==="
echo ""

SAMPLE_LIST="node_${SLURM_ARRAY_TASK_ID}_samples.txt"
printf '%s\n' "${NODE_SAMPLES[@]}" > "${SAMPLE_LIST}"

parallel -j 2 --progress --joblog "parallel_log_node${SLURM_ARRAY_TASK_ID}.txt" \
    count_sample :::: "${SAMPLE_LIST}"

rm -f "${SAMPLE_LIST}"

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=== Node ${SLURM_ARRAY_TASK_ID} Complete ==="

SUCCESS=0
ZERO_PATHOGEN=0
for bam_file in "${NODE_SAMPLES[@]}"; do
    sample_name=$(basename "${bam_file}" | sed 's/_Aligned.sortedByCoord.out.bam//')

    if [[ ${bam_file} == *"/soybean_samples/"* ]]; then
        output="${COUNTS_DIR}/soybean/${sample_name}_combined_counts.txt"
        path_summary="${COUNTS_DIR}/soybean/${sample_name}_pathogen_summary.txt"
    else
        output="${COUNTS_DIR}/sunflower/${sample_name}_combined_counts.txt"
        path_summary="${COUNTS_DIR}/sunflower/${sample_name}_pathogen_summary.txt"
    fi

    if [ -f "${output}" ]; then
        ((SUCCESS++))
        if [ -f "${path_summary}" ]; then
            assigned=$(awk '$1=="Assigned"{print $2}' "${path_summary}")
            if [ "${assigned:-0}" -eq 0 ]; then
                ((ZERO_PATHOGEN++))
            fi
        fi
    fi
done

printf "%-25s %d\n" "Assigned:" ${#NODE_SAMPLES[@]}
printf "%-25s %d\n" "Success:"  ${SUCCESS}
printf "%-25s %d\n" "Zero pathogen reads:" ${ZERO_PATHOGEN}
echo ""

echo "End time: $(date)"
ELAPSED=$SECONDS
echo "Runtime: $(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"

exit 0
