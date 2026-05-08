#!/bin/bash -l
# CRAM file renaming and CSV metadata generation script
# This script creates a CSV mapping and optionally renames CRAM files

#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --tmp=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wan00965@umn.edu
#SBATCH -e logs/cram_rename.err
#SBATCH -o logs/cram_rename.out
#SBATCH --job-name=cram_rename

# Configuration
SOURCE_DIR="/scratch.global/wan00965/Quick_test/data/seqs"
OUTPUT_CSV="/scratch.global/wan00965/Quick_test/analysis/sample_log.csv"
DRY_RUN=false  # IMPORTANT: Must be lowercase! Set to false to actually rename files
BACKUP_DIR="/scratch.global/wan00965/Quick_test/data/seqs/original_cram_backup"

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_CSV")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/rename_log_$(date +%Y%m%d_%H%M%S).txt"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to log messages
log_message() {
    echo -e "$1" | tee -a "$LOG_FILE"
}

# Function to extract metadata from filename
extract_metadata() {
    local filename=$1
    local isolate=$2
    local original_file=$3
    
    # Extract time point (24, 48, or 96) - will be 0 for NC files
    local time_point=$(echo "$filename" | grep -oE '^[0-9]+' | head -1)
    
    # Extract host (Gm or Ha) - check both filename and original file
    local host_code=""
    if [[ "$filename" =~ _Gm_ ]] || [[ "$original_file" =~ -Gm_ ]] || [[ "$original_file" =~ _Gm_ ]]; then
        host_code="Gm"
    elif [[ "$filename" =~ _Ha_ ]] || [[ "$original_file" =~ -Ha_ ]] || [[ "$original_file" =~ _Ha_ ]]; then
        host_code="Ha"
    fi
    
    # Convert host code to full name
    local host="unknown"
    case "$host_code" in
        "Gm") host="soybean" ;;
        "Ha") host="sunflower" ;;
    esac
    
    # Extract replicate number
    local replicate=""
    if [[ "$isolate" == "NC" ]]; then
        # For NC: extract number from filename pattern 0_X_...
        replicate=$(echo "$filename" | cut -d'_' -f2)
    else
        # For regular samples: extract second number (after timepoint)
        replicate=$(echo "$filename" | cut -d'_' -f2)
    fi
    
    # Determine aggressiveness and host specificity (combined field)
    local agg_host_spec="unknown"
    
    # Normalize isolate name (handle SSPotter vs SsPotter)
    local normalized_isolate=$(echo "$isolate" | sed 's/SSPotter/SsPotter/g')
    
    case "$normalized_isolate" in
        "WISS47"|"JS659")
            agg_host_spec="low_general"
            ;;
        "Xtra7"|"MNSS6")
            agg_host_spec="high_general"
            ;;
        "SsPotter")
            agg_host_spec="host_specific_soybean"
            ;;
        "BN172")
            agg_host_spec="host_specific_sunflower"
            ;;
        "MNSS4"|"SSMM2")
            agg_host_spec="low_general"
            ;;
        "NC")
            agg_host_spec="negative_control"
            time_point="0"
            ;;
        *)
            if [[ "$isolate" == *"Control"* ]] || [[ "$isolate" == *"Neg"* ]] || [[ "$isolate" == "NC" ]]; then
                agg_host_spec="negative_control"
                time_point="0"
            fi
            ;;
    esac
    
    echo "$host,$normalized_isolate,$time_point,$agg_host_spec,$replicate"
}

# Function to generate new filename
generate_new_filename() {
    local original=$1
    
    # Check if this is a negative control file
    if [[ "$original" =~ _NC_ ]]; then
        # Extract host (Gm or Ha)
        local host=$(echo "$original" | grep -oE '(Gm|Ha)' | head -1)
        # Extract NC number
        local nc_num=$(echo "$original" | grep -oE 'NC_[0-9]+' | grep -oE '[0-9]+')
        # Pattern: 424835-Gm_NC_1-Z0349-...cram → 0_1_Gm_NC.cram
        local new_name="0_${nc_num}_${host}_NC.cram"
    else
        # Regular sample
        # Pattern: {number}-{time}_{rep}_{host}_{isolate}-Z{num}-{barcode}.cram
        # Result: {time}_{rep}_{host}_{isolate}.cram
        local new_name=$(echo "$original" | sed -E 's/^[0-9]+-//; s/-Z[0-9]+-.*\.cram$/.cram/')
    fi
    
    echo "$new_name"
}

# Function to extract isolate name
extract_isolate() {
    local filename=$1
    
    # Check if this is a negative control
    if [[ "$filename" =~ _NC_ ]]; then
        echo "NC"
    else
        # Extract isolate name (last part before -Z)
        local isolate=$(echo "$filename" | grep -oP '(?<=_)[A-Za-z0-9]+(?=-Z)')
        # Normalize SSPotter to SsPotter
        echo "$isolate" | sed 's/SSPotter/SsPotter/g'
    fi
}

# Main function
main() {
    log_message "${GREEN}=========================================${NC}"
    log_message "${GREEN}CRAM File Renaming and Metadata Generator${NC}"
    log_message "${GREEN}=========================================${NC}"
    log_message "Start time: $(date)"
    log_message ""
    
    # Display dry run status prominently
    if [ "$DRY_RUN" = "true" ] || [ "$DRY_RUN" = "True" ] || [ "$DRY_RUN" = "TRUE" ]; then
        log_message "${YELLOW}╔════════════════════════════════════════╗${NC}"
        log_message "${YELLOW}║          DRY RUN MODE ACTIVE          ║${NC}"
        log_message "${YELLOW}║      NO FILES WILL BE RENAMED         ║${NC}"
        log_message "${YELLOW}║  Only CSV will be generated           ║${NC}"
        log_message "${YELLOW}╚════════════════════════════════════════╝${NC}"
        log_message ""
        log_message "To actually rename files: Set DRY_RUN=false in the script"
        log_message ""
        DRY_RUN="true"  # Normalize
    elif [ "$DRY_RUN" = "false" ] || [ "$DRY_RUN" = "False" ] || [ "$DRY_RUN" = "FALSE" ]; then
        log_message "${RED}╔════════════════════════════════════════╗${NC}"
        log_message "${RED}║          LIVE MODE ACTIVE             ║${NC}"
        log_message "${RED}║       FILES WILL BE RENAMED!          ║${NC}"
        log_message "${RED}╚════════════════════════════════════════╝${NC}"
        log_message ""
        DRY_RUN="false"  # Normalize
    else
        log_message "${RED}ERROR: DRY_RUN must be 'true' or 'false' (got: '$DRY_RUN')${NC}"
        exit 1
    fi
    
    # Verify source directory exists
    if [ ! -d "$SOURCE_DIR" ]; then
        log_message "${RED}ERROR: Source directory does not exist: $SOURCE_DIR${NC}"
        exit 1
    fi
    
    # Change to source directory
    cd "$SOURCE_DIR" || {
        log_message "${RED}ERROR: Cannot change to source directory${NC}"
        exit 1
    }
    
    log_message "${BLUE}Source directory: $SOURCE_DIR${NC}"
    log_message "${BLUE}Output CSV: $OUTPUT_CSV${NC}"
    log_message "${BLUE}Log file: $LOG_FILE${NC}"
    log_message ""
    
    # Check if running in dry run mode
    if [ "$DRY_RUN" = "false" ]; then
        # Create backup directory
        mkdir -p "$BACKUP_DIR"
        log_message "Backup directory: $BACKUP_DIR"
        log_message ""
    fi
    
    # Count total files
    total_files=$(ls *.cram 2>/dev/null | wc -l)
    
    if [ $total_files -eq 0 ]; then
        log_message "${RED}ERROR: No CRAM files found in $SOURCE_DIR${NC}"
        exit 1
    fi
    
    # Count regular samples vs controls
    regular_samples=$(ls *.cram 2>/dev/null | grep -v "_NC_" | wc -l)
    control_samples=$(ls *.cram 2>/dev/null | grep "_NC_" | wc -l)
    
    log_message "${GREEN}Found $total_files CRAM files to process${NC}"
    log_message "  Regular samples: $regular_samples"
    log_message "  Negative controls: $control_samples"
    log_message ""
    
    # Create CSV header matching Excel format
    echo "sample_id,file_name,Plant_host,Ss_isolate,time_point,aggressiveness_n_host_specificity,replicate" > "$OUTPUT_CSV"
    
    # Counter
    local count=0
    local error_count=0
    local skip_count=0
    
    # Process all CRAM files
    for cram_file in *.cram; do
        if [ ! -f "$cram_file" ]; then
            continue
        fi
        
        ((count++))
        
        # Generate new filename
        new_filename=$(generate_new_filename "$cram_file")
        
        # Extract isolate name
        isolate=$(extract_isolate "$cram_file")
        
        if [ -z "$isolate" ]; then
            log_message "${RED}[$count/$total_files] ERROR: Could not extract isolate from $cram_file${NC}"
            ((error_count++))
            continue
        fi
        
        # Extract metadata (pass original filename for host detection)
        metadata=$(extract_metadata "$new_filename" "$isolate" "$cram_file")
        
        # Create sample_id (same as new_filename without .cram)
        sample_id="${new_filename%.cram}"
        
        # Write to CSV in format: sample_id,file_name,Plant_host,Ss_isolate,time_point,aggressiveness_n_host_specificity,replicate
        echo "$sample_id,$cram_file,$metadata" >> "$OUTPUT_CSV"
        
        # Display progress (more concise for large batches)
        if [[ "$cram_file" =~ _NC_ ]]; then
            log_message "${BLUE}[$count/$total_files]${NC} $cram_file ${YELLOW}(Control)${NC}"
        else
            log_message "${GREEN}[$count/$total_files]${NC} $cram_file"
        fi
        log_message "  -> ${BLUE}$new_filename${NC}"
        
        # Rename file if not in dry run mode
        if [ "$DRY_RUN" = "false" ]; then
            if [ -f "$new_filename" ]; then
                log_message "${YELLOW}  WARNING: $new_filename already exists, skipping...${NC}"
                ((skip_count++))
            else
                # Copy to backup first
                if cp "$cram_file" "$BACKUP_DIR/"; then
                    # Rename
                    if mv "$cram_file" "$new_filename"; then
                        log_message "${GREEN}  ✓ Renamed${NC}"
                    else
                        log_message "${RED}  ✗ Failed to rename${NC}"
                        ((error_count++))
                    fi
                else
                    log_message "${RED}  ✗ Failed to backup${NC}"
                    ((error_count++))
                fi
            fi
        fi
        
        log_message ""
    done
    
    # Summary
    log_message ""
    log_message "${GREEN}=========================================${NC}"
    log_message "${GREEN}Summary${NC}"
    log_message "${GREEN}=========================================${NC}"
    log_message "End time: $(date)"
    log_message "Total files processed: $count / $total_files"
    log_message "  Regular samples: $regular_samples"
    log_message "  Negative controls: $control_samples"
    
    if [ "$DRY_RUN" = "false" ]; then
        log_message "Successfully renamed: $((count - error_count - skip_count))"
        log_message "Skipped (already exist): $skip_count"
    fi
    
    log_message "Errors: $error_count"
    log_message "CSV file created: $OUTPUT_CSV"
    log_message ""
    
    if [ "$DRY_RUN" = "true" ]; then
        log_message "${YELLOW}╔════════════════════════════════════════╗${NC}"
        log_message "${YELLOW}║     THIS WAS A DRY RUN                ║${NC}"
        log_message "${YELLOW}║     No files were renamed             ║${NC}"
        log_message "${YELLOW}╚════════════════════════════════════════╝${NC}"
        log_message ""
        log_message "Next steps:"
        log_message "1. Review the CSV file: $OUTPUT_CSV"
        log_message "2. If everything looks correct, edit script:"
        log_message "   Change: DRY_RUN=true  to  DRY_RUN=false"
        log_message "3. Run again: bash $(basename $0)"
        log_message "   Or: sbatch $(basename $0)"
    else
        log_message "${GREEN}Files have been renamed!${NC}"
        log_message "Original files backed up to: $BACKUP_DIR"
    fi
    
    log_message "${GREEN}=========================================${NC}"
    
    # Display first few lines of CSV
    log_message ""
    log_message "${BLUE}First 10 rows of CSV:${NC}"
    head -11 "$OUTPUT_CSV" | tee -a "$LOG_FILE"
    
    log_message ""
    log_message "${BLUE}Last 10 rows of CSV (including controls):${NC}"
    tail -10 "$OUTPUT_CSV" | tee -a "$LOG_FILE"
}

# Run main function
main

exit 0