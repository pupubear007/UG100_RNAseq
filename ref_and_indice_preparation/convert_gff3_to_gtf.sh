#!/bin/bash
#SBATCH --job-name=gff3_to_gtf
#SBATCH --time=30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH -o logs/gff3_to_gtf_%j.out
#SBATCH -e logs/gff3_to_gtf_%j.err


# Convert GFF3 to GTF format for featureCounts compatibility

set -e

echo "=== Converting GFF3 to GTF Format ==="
echo ""

BASE_DIR="/scratch.global/wan00965/Quick_test"

# Install gffread if not available
module load conda
source activate pu_RNAseq

if ! command -v gffread &> /dev/null; then
    echo "Installing gffread..."
    conda install -c bioconda gffread -y
fi

# Convert soybean GFF3
SOYBEAN_GFF="${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gff3"
SOYBEAN_GTF="${BASE_DIR}/jgi_processed/combined/soybean_ss/Gm_Ss_genes.gtf"

if [ -f "${SOYBEAN_GFF}" ]; then
    echo "Converting soybean GFF3 to GTF..."
    gffread "${SOYBEAN_GFF}" -T -o "${SOYBEAN_GTF}"
    echo "  ✓ Created: ${SOYBEAN_GTF}"
    echo "  Size: $(ls -lh ${SOYBEAN_GTF} | awk '{print $5}')"
else
    echo "  ✗ Soybean GFF3 not found"
fi

echo ""

# Convert sunflower GFF3
SUNFLOWER_GFF="${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gff3"
SUNFLOWER_GTF="${BASE_DIR}/jgi_processed/combined/sunflower_ss/Ha_Ss_genes.gtf"

if [ -f "${SUNFLOWER_GFF}" ]; then
    echo "Converting sunflower GFF3 to GTF..."
    gffread "${SUNFLOWER_GFF}" -T -o "${SUNFLOWER_GTF}"
    echo "  ✓ Created: ${SUNFLOWER_GTF}"
    echo "  Size: $(ls -lh ${SUNFLOWER_GTF} | awk '{print $5}')"
else
    echo "  ✗ Sunflower GFF3 not found"
fi

echo ""
echo "✅ Conversion complete!"
echo ""
echo "GTF files created - these will work better with featureCounts."
echo "The alignment script will automatically use .gtf if available, otherwise .gff3"
