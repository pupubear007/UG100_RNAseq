#!/bin/bash
#SBATCH --job-name=collect_counts_salmon
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/salmon_count_only_%j.out
#SBATCH -e logs/salmon_count_only_%j.err


# ============================================================================
# COLLECT GENE COUNTS FROM SALMON ONLY
# Extracts gene-level counts from Salmon quant.genes.sf files
# ============================================================================

set -e

echo "=== Collecting Salmon Gene Counts ==="
echo "Start time: $(date)"
echo ""

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR="/scratch.global/wan00965/Quick_test"
RESULTS_DIR="${BASE_DIR}/analysis/final_results"
SALMON_DIR="${BASE_DIR}/analysis/salmon_quant_combined_ref"

# Create results directory
mkdir -p "${RESULTS_DIR}"

# ============================================================================
# Function to collect Salmon gene counts
# ============================================================================
collect_salmon_counts() {
    local species=$1
    local output_file="${RESULTS_DIR}/gene_counts_${species}_SALMON.txt"
    local host_only="${RESULTS_DIR}/gene_counts_${species}_SALMON_HOST_only.txt"
    local pathogen_only="${RESULTS_DIR}/gene_counts_${species}_SALMON_PATHOGEN_only.txt"
    local stats_file="${RESULTS_DIR}/pathogen_load_${species}_SALMON.txt"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Salmon: Collecting ${species} gene counts"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Find all Salmon quantification directories
    local salmon_dirs=($(find "${SALMON_DIR}/${species}_samples" -mindepth 1 -maxdepth 1 -type d | sort))
    local num_samples=${#salmon_dirs[@]}
    
    if [ ${num_samples} -eq 0 ]; then
        echo "  ✗ No Salmon quantification directories found for ${species}"
        echo "  Location checked: ${SALMON_DIR}/${species}_samples"
        return 1
    fi
    
    echo "  Found ${num_samples} samples"
    
    # Check if gene-level counts exist (quant.genes.sf from --geneMap option)
    local first_gene_file="${salmon_dirs[0]}/quant.genes.sf"
    
    if [ ! -f "${first_gene_file}" ]; then
        echo "  ✗ Gene-level quantification not found (quant.genes.sf)"
        echo "    Salmon was likely run without --geneMap parameter"
        echo "    Only transcript-level counts available (quant.sf)"
        echo ""
        echo "  To fix: Re-run alignment with --geneMap option"
        return 1
    fi
    
    echo "  ✓ Gene-level quantification files found"
    
    # Build header
    echo -ne "gene_id" > "${output_file}"
    
    for salmon_dir in "${salmon_dirs[@]}"; do
        sample_name=$(basename "${salmon_dir}")
        echo -ne "\t${sample_name}" >> "${output_file}"
    done
    echo "" >> "${output_file}"
    
    # Extract gene IDs from first file
    tail -n +2 "${first_gene_file}" | cut -f1 > "${output_file}.genes"
    
    # Extract NumReads (column 5) from all samples
    for salmon_dir in "${salmon_dirs[@]}"; do
        gene_file="${salmon_dir}/quant.genes.sf"
        tail -n +2 "${gene_file}" | cut -f5 > "${output_file}.tmp.$(basename ${salmon_dir})"
    done
    
    # Combine all count columns
    paste "${output_file}.genes" ${output_file}.tmp.* >> "${output_file}"
    
    # Cleanup temporary files
    rm "${output_file}.genes" ${output_file}.tmp.*
    
    # Separate host and pathogen genes
    if [[ ${species} == "soybean" ]]; then
        head -1 "${output_file}" > "${host_only}"
        grep "^Gm_" "${output_file}" >> "${host_only}"
        
        head -1 "${output_file}" > "${pathogen_only}"
        grep "^Ss_" "${output_file}" >> "${pathogen_only}"
    else
        head -1 "${output_file}" > "${host_only}"
        grep "^Ha_" "${output_file}" >> "${host_only}"
        
        head -1 "${output_file}" > "${pathogen_only}"
        grep "^Ss_" "${output_file}" >> "${pathogen_only}"
    fi
    
    # Calculate pathogen load per sample
    echo -e "sample\thost_reads\tpathogen_reads\tpathogen_load_percent" > "${stats_file}"
    
    num_cols=$(head -1 "${output_file}" | awk '{print NF}')
    
    for ((col=2; col<=num_cols; col++)); do
        sample_name=$(head -1 "${output_file}" | cut -f${col})
        
        if [[ ${species} == "soybean" ]]; then
            host_reads=$(grep "^Gm_" "${output_file}" | cut -f${col} | awk '{sum+=$1} END {printf "%.0f", sum}')
            pathogen_reads=$(grep "^Ss_" "${output_file}" | cut -f${col} | awk '{sum+=$1} END {printf "%.0f", sum}')
        else
            host_reads=$(grep "^Ha_" "${output_file}" | cut -f${col} | awk '{sum+=$1} END {printf "%.0f", sum}')
            pathogen_reads=$(grep "^Ss_" "${output_file}" | cut -f${col} | awk '{sum+=$1} END {printf "%.0f", sum}')
        fi
        
        total_reads=$((host_reads + pathogen_reads))
        if [ ${total_reads} -gt 0 ]; then
            pathogen_load=$(echo "scale=4; ${pathogen_reads} / ${total_reads} * 100" | bc)
        else
            pathogen_load=0
        fi
        
        echo -e "${sample_name}\t${host_reads}\t${pathogen_reads}\t${pathogen_load}" >> "${stats_file}"
    done
    
    local total_genes=$(tail -n +2 "${output_file}" | wc -l)
    local host_genes=$(tail -n +2 "${host_only}" | wc -l)
    local pathogen_genes=$(tail -n +2 "${pathogen_only}" | wc -l)
    
    echo "  ✓ Total genes: ${total_genes}"
    echo "  ✓ Host genes: ${host_genes}"
    echo "  ✓ Pathogen genes: ${pathogen_genes}"
    echo "  ✓ Output: ${output_file}"
    echo ""
}

# ============================================================================
# Process both species
# ============================================================================

echo "=== Processing Soybean Samples ==="
echo ""
collect_salmon_counts "soybean"

echo "=== Processing Sunflower Samples ==="
echo ""
collect_salmon_counts "sunflower"

# ============================================================================
# Generate summary report
# ============================================================================
SUMMARY="${RESULTS_DIR}/SALMON_COLLECTION_SUMMARY.txt"

echo "=== Salmon Gene Count Collection Summary ===" > "${SUMMARY}"
echo "Generated: $(date)" >> "${SUMMARY}"
echo "" >> "${SUMMARY}"

echo "Output Location: ${RESULTS_DIR}/" >> "${SUMMARY}"
echo "" >> "${SUMMARY}"

echo "Files Created:" >> "${SUMMARY}"
ls -lh "${RESULTS_DIR}/"*SALMON* 2>/dev/null | tail -n +2 | awk '{printf "  %-50s %8s\n", $9, $5}' >> "${SUMMARY}"

echo "" >> "${SUMMARY}"
echo "✅ Salmon gene count collection complete!" >> "${SUMMARY}"

cat "${SUMMARY}"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ SALMON GENE COUNTS COLLECTED SUCCESSFULLY!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 Results available in: ${RESULTS_DIR}/"
echo ""
echo "Files created:"
echo "  • gene_counts_soybean_SALMON.txt"
echo "  • gene_counts_sunflower_SALMON.txt"
echo "  • gene_counts_*_SALMON_HOST_only.txt"
echo "  • gene_counts_*_SALMON_PATHOGEN_only.txt"
echo "  • pathogen_load_soybean_SALMON.txt"
echo "  • pathogen_load_sunflower_SALMON.txt"
echo ""
echo "Next step:"
echo "  Use gene_counts_*_SALMON.txt as input for DESeq2 analysis"
echo ""
echo "End time: $(date)"

exit 0
