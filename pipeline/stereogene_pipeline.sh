#!/bin/bash

# Define directories
alignments_dir="./alignments"
p_value_dir="./p-value_peaks"
fold_dir="./fold_change_over_control"
alignments_dir="./alignments"

# Create list of files in each directory and save it to list.lst
for dir in "$alignments_dir" "$control_dir" "$p_value_dir" "$fold_dir" "$alignments_dir"; do
    echo "Processing directory: $dir"
    
    # Create list of files
    ls -1 "$dir" > "$dir/list.lst"
    
    # Run commands for each directory
    cd "$dir" || exit
    
    # Run StereoGene, Confounder, Projector, and StereoGene again
    ./StereoGene cfg=../stereogene.cfg list.lst
    ./Confounder cfg=../stereogene.cfg list.lst
    ./Projector cfg=../confounder.cfg list.lst
    ./StereoGene cfg=../stereogene.cfg list.lst
    
    cd - # Return to the original directory
done

echo "Processing completed!"
