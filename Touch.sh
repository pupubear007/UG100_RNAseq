#!/bin/bash
# touch_progress.sh

TARGET="${1:-/scratch.global/wan00965/Quick_test/}"
COUNT=0
TOTAL=$(find "$TARGET" -type f | wc -l)
echo "Found $TOTAL files. Starting..."

START=$(date +%s)

find "$TARGET" -type f -print0 | while IFS= read -r -d '' file; do
    touch -am "$file"
    COUNT=$((COUNT + 1))
    if (( COUNT % 100 == 0 )); then
        ELAPSED=$(( $(date +%s) - START ))
        printf "\rProcessed: %d / %d files (%ds elapsed)" "$COUNT" "$TOTAL" "$ELAPSED"
    fi
done
echo -e "\nDone."