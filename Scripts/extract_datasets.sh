#!/bin/bash

# Source directory (read-only)
SOURCE_DIR="/local/ofir/hadasra/create_sc_object_workflow/temp_objects"

# Your personal directory to save the output lists
OUTPUT_DIR="/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/dataset_per_cancer_lists"
# Create your output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List of cancer types to extract
CANCER_TYPES=("Breast" "Pancreas" "Lung" "Ovarian" "Kidney" "Skin" "Sarcoma" "Neuroendocrine" "Colorectal" "Brain" "Liver" "Prostate" "Head_and_Neck")
# Create a temporary file to track filenames for subtraction logic
MATCHED_FILENAMES=$(mktemp)

echo "Processing files from: $SOURCE_DIR"
echo "----------------------------------"

for TYPE in "${CANCER_TYPES[@]}"; do
    # 1. Save FULL PATHS to the type-specific file
    find "$SOURCE_DIR" -maxdepth 1 -iname "*${TYPE}*.RDS" > "$OUTPUT_DIR/${TYPE}_files.txt"

    # 2. Track just the FILENAMES (without paths) for the "unknown" subtraction logic later
    if [ -s "$OUTPUT_DIR/${TYPE}_files.txt" ]; then
        find "$SOURCE_DIR" -maxdepth 1 -iname "*${TYPE}*.RDS" -printf "%f\n" >> "$MATCHED_FILENAMES"
        
        COUNT=$(wc -l < "$OUTPUT_DIR/${TYPE}_files.txt")
        echo "Found $COUNT .RDS files for $TYPE (Full paths saved)."
    else
        rm "$OUTPUT_DIR/${TYPE}_files.txt"
    fi
done

# Check for unmatched files - using -printf "%f\n" to keep only names as requested
echo "Checking for unmatched files..."
find "$SOURCE_DIR" -maxdepth 1 -iname "*.RDS" -printf "%f\n" | grep -vFf "$MATCHED_FILENAMES" > "$OUTPUT_DIR/unknown_files.txt"

# Final cleanup for the unknown list
if [ ! -s "$OUTPUT_DIR/unknown_files.txt" ]; then
    rm "$OUTPUT_DIR/unknown_files.txt"
    echo "No unknown files found."
else
    UNKNOWN_COUNT=$(wc -l < "$OUTPUT_DIR/unknown_files.txt")
    echo "Found $UNKNOWN_COUNT unknown .RDS files (Filenames only saved to unknown_files.txt)."
fi

rm "$MATCHED_FILENAMES"
echo "----------------------------------"
echo "Done."