#!/bin/bash
# create_count_matrix_FIXED.sh - Properly merge counts with gene ID matching

BASE_DIR="/scratch.global/wan00965/Quick_test"
COUNTS_DIR="${BASE_DIR}/analysis/dual_counts_v5"
OUTPUT_DIR="${BASE_DIR}/analysis/count_matrices_fixed"

mkdir -p "${OUTPUT_DIR}"

echo "=== Creating Count Matrices (FIXED VERSION) ==="
echo ""

# ============================================================================
# Function to create matrix with proper gene matching
# ============================================================================
create_matrix() {
    local host_type=$1
    local input_dir="${COUNTS_DIR}/${host_type}"
    local output_file="${OUTPUT_DIR}/${host_type}_count_matrix.txt"
    local temp_dir="${OUTPUT_DIR}/temp_${host_type}"
    
    mkdir -p "${temp_dir}"
    
    echo "Processing ${host_type}..."
    
    # Get list of count files (sorted)
    COUNT_FILES=($(ls -1 "${input_dir}"/*_combined_counts.txt 2>/dev/null | sort))
    
    if [ ${#COUNT_FILES[@]} -eq 0 ]; then
        echo "  ERROR: No count files found in ${input_dir}"
        return 1
    fi
    
    echo "  Found ${#COUNT_FILES[@]} samples"
    
    # ========================================================================
    # Step 1: Extract and sort all individual files by gene ID
    # ========================================================================
    echo "  Sorting individual count files..."
    for i in "${!COUNT_FILES[@]}"; do
        count_file="${COUNT_FILES[$i]}"
        sample_name=$(basename "${count_file}" | sed 's/_combined_counts.txt//')
        
        # Extract gene ID and count, skip header lines, remove quotes, sort by gene ID
        tail -n +3 "${count_file}" | \
            sed 's/"//g' | \
            awk '{print $1"\t"$2}' | \
            sort -k1,1 > "${temp_dir}/${sample_name}_sorted.txt"
        
        printf "\r  Progress: %d/%d samples" $((i+1)) ${#COUNT_FILES[@]}
    done
    echo ""
    
    # ========================================================================
    # Step 2: Get complete list of unique genes (union of all genes)
    # ========================================================================
    echo "  Building gene list..."
    cat "${temp_dir}"/*_sorted.txt | \
        cut -f1 | \
        sort -u > "${temp_dir}/all_genes.txt"
    
    TOTAL_GENES=$(wc -l < "${temp_dir}/all_genes.txt")
    echo "  ${TOTAL_GENES} unique genes detected"
    
    # ========================================================================
    # Step 3: Create header
    # ========================================================================
    echo -ne "GeneID" > "${output_file}"
    for count_file in "${COUNT_FILES[@]}"; do
        sample_name=$(basename "${count_file}" | sed 's/_combined_counts.txt//')
        echo -ne "\t${sample_name}" >> "${output_file}"
    done
    echo "" >> "${output_file}"
    
    # ========================================================================
    # Step 4: Merge counts using join (proper gene ID matching)
    # ========================================================================
    echo "  Merging counts with gene ID matching..."
    
    # Start with gene list
    cp "${temp_dir}/all_genes.txt" "${temp_dir}/matrix.txt"
    
    # Join each sample's counts
    for i in "${!COUNT_FILES[@]}"; do
        count_file="${COUNT_FILES[$i]}"
        sample_name=$(basename "${count_file}" | sed 's/_combined_counts.txt//')
        
        # Join on gene ID, fill missing with 0
        join -a 1 -1 1 -2 1 -o auto -e "0" -t $'\t' \
            "${temp_dir}/matrix.txt" \
            "${temp_dir}/${sample_name}_sorted.txt" \
            > "${temp_dir}/matrix_temp.txt"
        
        mv "${temp_dir}/matrix_temp.txt" "${temp_dir}/matrix.txt"
        
        printf "\r  Progress: %d/%d samples" $((i+1)) ${#COUNT_FILES[@]}
    done
    echo ""
    
    # ========================================================================
    # Step 5: Output final matrix
    # ========================================================================
    echo "  Writing final matrix..."
    cat "${temp_dir}/matrix.txt" >> "${output_file}"
    
    # Cleanup
    rm -rf "${temp_dir}"
    
    # Summary
    ROWS=$(wc -l < "${output_file}")
    COLS=$(head -1 "${output_file}" | awk '{print NF}')
    
    echo "  ✓ Matrix created: $((ROWS-1)) genes × $((COLS-1)) samples"
    echo "  Output: ${output_file}"
    echo ""
}

# ============================================================================
# Create matrices
# ============================================================================

create_matrix "soybean"
create_matrix "sunflower"

echo "=== Summary ==="
echo ""
ls -lh "${OUTPUT_DIR}"/*.txt
echo ""
echo "Done! Count matrices ready for analysis."
echo ""
echo "IMPORTANT: These fixed matrices are in:"
echo "  ${OUTPUT_DIR}/"
echo ""
echo "Compare with original matrices to verify the fix!"
