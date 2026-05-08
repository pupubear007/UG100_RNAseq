#!/bin/bash

#SBATCH --job-name=regen_sample_log
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/regenerate_sample_log_%j.out
#SBATCH -e logs/regenerate_sample_log_%j.err

# Script to regenerate sample_log.csv from renamed CRAM files
# Usage: sbatch regenerate_sample_log.sh

# Set paths
CRAM_DIR="/scratch.global/wan00965/Quick_test/data/seqs"
OUTPUT_DIR="/scratch.global/wan00965/Quick_test/analysis"
OUTPUT_CSV="${OUTPUT_DIR}/sample_log.csv"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "=== Regenerating Sample Log CSV ==="
echo "CRAM directory: ${CRAM_DIR}"
echo "Output CSV: ${OUTPUT_CSV}"
echo ""

# Count CRAM files
CRAM_COUNT=$(find "${CRAM_DIR}" -maxdepth 1 -name "*.cram" | wc -l)
echo "Found ${CRAM_COUNT} CRAM files"
echo ""

# Create CSV header
echo "sample_id,file_name,Plant_host,Ss_isolate,time_point,aggressiveness_n_host_specificity,replicate" > "${OUTPUT_CSV}"

# Function to determine aggressiveness category
get_aggressiveness_category() {
    local isolate=$1
    
    case "${isolate}" in
        "WISS47"|"JS659")
            echo "low_general"
            ;;
        "Xtra7"|"MNSS6")
            echo "high_general"
            ;;
        "SsPotter"|"SSPotter")
            echo "host_specific_soybean"
            ;;
        "BN172")
            echo "host_specific_sunflower"
            ;;
        "NC")
            echo "negative_control"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

# Function to get full host name
get_host_name() {
    local host_code=$1
    
    case "${host_code}" in
        "Gm")
            echo "soybean"
            ;;
        "Ha")
            echo "sunflower"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

# Process each CRAM file
COUNTER=0
for cram_file in "${CRAM_DIR}"/*.cram; do
    # Skip if no files found
    [[ -e "$cram_file" ]] || continue
    
    # Get base filename without path and extension
    basename=$(basename "${cram_file}" .cram)
    
    # Parse filename: {time}_{rep}_{host}_{isolate}
    # Example: 24_1_Gm_WISS47 or 0_1_Gm_NC
    IFS='_' read -ra PARTS <<< "${basename}"
    
    TIME="${PARTS[0]}"
    REP="${PARTS[1]}"
    HOST_CODE="${PARTS[2]}"
    ISOLATE="${PARTS[3]}"
    
    # Get full host name
    HOST_NAME=$(get_host_name "${HOST_CODE}")
    
    # Get aggressiveness category
    AGGRESSIVENESS=$(get_aggressiveness_category "${ISOLATE}")
    
    # Create sample_id
    SAMPLE_ID="${basename}"
    
    # Write to CSV
    echo "${SAMPLE_ID},${basename}.cram,${HOST_NAME},${ISOLATE},${TIME},${AGGRESSIVENESS},${REP}" >> "${OUTPUT_CSV}"
    
    ((COUNTER++))
done

echo "=== Processing Complete ==="
echo "Processed ${COUNTER} CRAM files"
echo "CSV saved to: ${OUTPUT_CSV}"
echo ""

# Show first 10 lines of the CSV
echo "=== First 10 lines of sample_log.csv ==="
head -n 10 "${OUTPUT_CSV}"
echo ""

# Summary statistics
echo "=== Summary Statistics ==="
echo "Total samples: $((COUNTER))"
echo ""
echo "By Plant Host:"
tail -n +2 "${OUTPUT_CSV}" | cut -d',' -f3 | sort | uniq -c
echo ""
echo "By Isolate:"
tail -n +2 "${OUTPUT_CSV}" | cut -d',' -f4 | sort | uniq -c
echo ""
echo "By Time Point:"
tail -n +2 "${OUTPUT_CSV}" | cut -d',' -f5 | sort -n | uniq -c
echo ""
echo "By Aggressiveness Category:"
tail -n +2 "${OUTPUT_CSV}" | cut -d',' -f6 | sort | uniq -c
echo ""

echo "✅ Sample log regeneration complete!"