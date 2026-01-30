#!/bin/bash
# Directory containing the results
results_dir="/burg/cgl/users/mv2850/runs/250903_mv_jurkat_adept/mixcr/results"
export_dir="$results_dir/export"

# Create the export directory if it doesn't exist
mkdir -p "$export_dir"

# Find and copy all files with ".clones" in the filename
find "$results_dir" -type f | grep "\.clones" | while read -r file; do
    cp "$file" "$export_dir"
done
echo "All .clones files have been copied to $export_dir"
