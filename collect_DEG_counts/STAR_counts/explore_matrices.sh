#!/bin/bash
# explore_matrices.sh - Quick exploration of your count matrices

BASE_DIR="/scratch.global/wan00965/Quick_test"
MATRIX_DIR="${BASE_DIR}/analysis/count_matrices_fixed"

echo "=== Count Matrix Exploration ==="
echo ""

cd "${MATRIX_DIR}" || exit 1

# ============================================================================
# Soybean Matrix
# ============================================================================
echo "📊 SOYBEAN MATRIX"
echo "================="
echo ""

SOYBEAN="soybean_count_matrix.txt"

echo "Dimensions:"
ROWS=$(wc -l < "${SOYBEAN}")
COLS=$(head -1 "${SOYBEAN}" | awk '{print NF}')
echo "  Genes: $((ROWS - 1))"
echo "  Samples: $((COLS - 1))"
echo ""

echo "Sample names (first 5):"
head -1 "${SOYBEAN}" | cut -f2-6 | tr '\t' '\n' | head -5
echo "  ..."
echo ""

echo "Gene types:"
HOST_GENES=$(grep -c "^Glyma" "${SOYBEAN}")
PATH_GENES=$(grep -c "^SS1G" "${SOYBEAN}")
echo "  Host (Glyma.*): ${HOST_GENES}"
echo "  Pathogen (SS1G_*): ${PATH_GENES}"
echo ""

echo "Library sizes (total counts per sample):"
awk 'NR>1 {for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) print sum[i]}' "${SOYBEAN}" | \
    awk '{sum+=$1; n++} END {printf "  Min: %d\n  Max: %d\n  Mean: %d\n", min, max, sum/n}' \
        min="$(awk 'NR>1 {for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) print sum[i]}' "${SOYBEAN}" | sort -n | head -1)" \
        max="$(awk 'NR>1 {for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) print sum[i]}' "${SOYBEAN}" | sort -n | tail -1)"
echo ""

echo "Pathogen gene counts (should see higher counts in infected samples):"
grep "^SS1G" "${SOYBEAN}" | awk 'NR==1 {for(i=2;i<=NF;i++) sum[i]=$i} {for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) printf "%d ", sum[i]; print ""}' | \
    tr ' ' '\n' | grep -v "^$" | sort -n | awk '{arr[NR]=$1} END {
        printf "  Min: %d\n  Median: %d\n  Max: %d\n", arr[1], arr[int(NR/2)], arr[NR]
    }'
echo ""

# ============================================================================
# Sunflower Matrix
# ============================================================================
echo "📊 SUNFLOWER MATRIX"
echo "==================="
echo ""

SUNFLOWER="sunflower_count_matrix.txt"

echo "Dimensions:"
ROWS=$(wc -l < "${SUNFLOWER}")
COLS=$(head -1 "${SUNFLOWER}" | awk '{print NF}')
echo "  Genes: $((ROWS - 1))"
echo "  Samples: $((COLS - 1))"
echo ""

echo "Sample names (first 5):"
head -1 "${SUNFLOWER}" | cut -f2-6 | tr '\t' '\n' | head -5
echo "  ..."
echo ""

echo "Gene types:"
HOST_GENES=$(grep -c "^HanXRQ" "${SUNFLOWER}")
PATH_GENES=$(grep -c "^SS1G" "${SUNFLOWER}")
echo "  Host (HanXRQ*): ${HOST_GENES}"
echo "  Pathogen (SS1G_*): ${PATH_GENES}"
echo ""

echo "Library sizes (total counts per sample):"
awk 'NR>1 {for(i=2;i<=NF;i++) sum[i]+=$i} END {for(i=2;i<=NF;i++) print sum[i]}' "${SUNFLOWER}" | \
    awk '{sum+=$1; n++; if(NR==1 || $1<min) min=$1; if(NR==1 || $1>max) max=$1} END {printf "  Min: %d\n  Max: %d\n  Mean: %d\n", min, max, sum/n}'
echo ""

# ============================================================================
# Split matrices by host/pathogen
# ============================================================================
echo "=== Creating Split Matrices ==="
echo ""

echo "Creating soybean_host_only.txt..."
head -1 "${SOYBEAN}" > soybean_host_only.txt
grep "^Glyma" "${SOYBEAN}" >> soybean_host_only.txt
echo "  ✓ $(wc -l < soybean_host_only.txt) genes"

echo "Creating pathogen_soybean.txt..."
head -1 "${SOYBEAN}" > pathogen_soybean.txt
grep "^SS1G" "${SOYBEAN}" >> pathogen_soybean.txt
echo "  ✓ $(wc -l < pathogen_soybean.txt) genes"

echo "Creating sunflower_host_only.txt..."
head -1 "${SUNFLOWER}" > sunflower_host_only.txt
grep "^HanXRQ" "${SUNFLOWER}" >> sunflower_host_only.txt
echo "  ✓ $(wc -l < sunflower_host_only.txt) genes"

echo "Creating pathogen_sunflower.txt..."
head -1 "${SUNFLOWER}" > pathogen_sunflower.txt
grep "^SS1G" "${SUNFLOWER}" >> pathogen_sunflower.txt
echo "  ✓ $(wc -l < pathogen_sunflower.txt) genes"

echo ""
echo "=== Files Created ==="
ls -lh *.txt
echo ""
echo "Done! Ready for analysis in R/Python."
