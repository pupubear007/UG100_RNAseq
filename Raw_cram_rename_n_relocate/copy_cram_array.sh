#!/bin/bash -l

#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -e logs/movefile_%A_%a.err
#SBATCH -o logs/movefile_%A_%a.out
#SBATCH --job-name=copy_cram_array
#SBATCH --array=0-119  # Adjust based on number of files (e.g., 0-119 for 120 files)

num_files=$(ls /projects/standard/mmccaghe/shared/SS_Isolate_aggressiveness_RNAseq/*.cram | wc -l)
echo "Found $num_files files"

# Source and destination directories
SOURCE_DIR="/projects/standard/mmccaghe/shared/SS_Isolate_aggressiveness_RNAseq"
DEST_DIR="/scratch.global/wan00965/Quick_test/data/seqs"

# Create destination directory
mkdir -p "$DEST_DIR"

# Get list of all CRAM files into an array
mapfile -t CRAM_FILES < <(ls "$SOURCE_DIR"/*.cram 2>/dev/null)

# Get the file for this array task
FILE_TO_COPY="${CRAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$FILE_TO_COPY" ]; then
    echo "No file assigned to task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

filename=$(basename "$FILE_TO_COPY")

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Copying: $filename"
echo "From: $SOURCE_DIR"
echo "To: $DEST_DIR"

# Copy the file
if cp "$FILE_TO_COPY" "$DEST_DIR/"; then
    echo "✓ Successfully copied: $filename"
    exit 0
else
    echo "✗ Failed to copy: $filename"
    exit 1
fi