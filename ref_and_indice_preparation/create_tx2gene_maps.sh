#!/bin/bash
#SBATCH --job-name=create_tx2gene
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -o logs/tx2gene_maps_%j.out
#SBATCH -e logs/tx2gene_maps_%j.err

# ============================================================================
# Create Transcript-to-Gene Mapping Files for Salmon
# ============================================================================

set -e

BASE_DIR="/scratch.global/wan00965/Quick_test"
SOYBEAN_GTF="${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gff3"
SUNFLOWER_GTF="${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gff3"

OUTPUT_DIR="${BASE_DIR}/jgi_processed/tx2gene"
mkdir -p "${OUTPUT_DIR}"

echo "=== Creating Transcript-to-Gene Mapping Files ==="
echo "Start time: $(date)"
echo ""

# ============================================================================
# Function to create tx2gene from GFF3
# ============================================================================
create_tx2gene() {
    local gtf_file=$1
    local output_file=$2
    local species=$3
    
    echo "Processing: ${species}"
    echo "  Input:  ${gtf_file}"
    echo "  Output: ${output_file}"
    
    # Extract transcript-to-gene mapping from GFF3
    # Format: transcript_id<TAB>gene_id
    awk -F'\t' '
        $3 == "mRNA" || $3 == "transcript" {
            # Extract ID and Parent from attributes (column 9)
            split($9, attrs, ";")
            transcript_id = ""
            gene_id = ""
            
            for (i in attrs) {
                if (attrs[i] ~ /^ID=/) {
                    gsub(/^ID=/, "", attrs[i])
                    transcript_id = attrs[i]
                }
                if (attrs[i] ~ /Parent=/) {
                    gsub(/Parent=/, "", attrs[i])
                    gene_id = attrs[i]
                }
            }
            
            if (transcript_id != "" && gene_id != "") {
                print transcript_id "\t" gene_id
            }
        }
    ' "${gtf_file}" > "${output_file}"
    
    local num_mappings=$(wc -l < "${output_file}")
    echo "  ✓ Created ${num_mappings} transcript-to-gene mappings"
    echo ""
}

# ============================================================================
# Create mapping files
# ============================================================================

if [ -f "${SOYBEAN_GTF}" ]; then
    create_tx2gene "${SOYBEAN_GTF}" "${OUTPUT_DIR}/soybean_ss_tx2gene.tsv" "Soybean + Ss"
else
    echo "✗ ERROR: Soybean GTF not found: ${SOYBEAN_GTF}"
    exit 1
fi

if [ -f "${SUNFLOWER_GTF}" ]; then
    create_tx2gene "${SUNFLOWER_GTF}" "${OUTPUT_DIR}/sunflower_ss_tx2gene.tsv" "Sunflower + Ss"
else
    echo "✗ ERROR: Sunflower GTF not found: ${SUNFLOWER_GTF}"
    exit 1
fi

# ============================================================================
# Verify outputs
# ============================================================================
echo "=== Verification ==="
echo ""

for species in soybean sunflower; do
    tx2gene="${OUTPUT_DIR}/${species}_ss_tx2gene.tsv"
    
    if [ -f "${tx2gene}" ]; then
        size=$(ls -lh "${tx2gene}" | awk '{print $5}')
        lines=$(wc -l < "${tx2gene}")
        
        echo "${species}_ss_tx2gene.tsv:"
        echo "  Size:  ${size}"
        echo "  Lines: ${lines}"
        echo "  Preview:"
        head -3 "${tx2gene}" | sed 's/^/    /'
        echo ""
    fi
done

echo "✅ Transcript-to-gene mapping files created!"
echo "Location: ${OUTPUT_DIR}/"
echo ""
echo "End time: $(date)"

exit 0
