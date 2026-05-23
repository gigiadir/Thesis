#!/bin/bash

# Define the target directory
TARGET_DIR="/gpfs0/bgu-ofircohen/users/gigiadir/CCC-scDiffCom/results/split-by-rank-genes"

# Navigate to the directory
cd "$TARGET_DIR" || { echo "Directory not found"; exit 1; }

# Define the datasets to group by
DATASETS=("HNSCC.Atlas" "Choi_HNSC" "Kurten_HNSC" "Puram_HNSC" "Bill_HNSC")
# Loop through each dataset name
for DATASET in "${DATASETS[@]}"; do
    # Create the directory if it doesn't exist
    mkdir -p "$DATASET"
    
    # Move files that contain the dataset string in their name
    # The * wildcard matches the prefix and suffix around the dataset name
    mv *"$DATASET"* "$DATASET"/ 2>/dev/null
    
    echo "Grouped files for: $DATASET"
done

echo "Done organizing your Thesis output files."