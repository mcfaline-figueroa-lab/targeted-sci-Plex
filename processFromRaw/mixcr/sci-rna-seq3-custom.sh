#!/bin/bash
#
#replace ACCOUNT with your account name before submitting.
#
#SBATCH --account=cgl            # Replace ACCOUNT with your group account name
#SBATCH --job-name=mixcr     # The job name
#SBATCH -c 16                     # The number of cores to request
#SBATCH --time=0-6:00            # The time the job will take to run in D-HH:MM
#SBATCH --mem-per-cpu=6G                # Memory limit

# Load necessary modules (if needed)
#conda activate /burg/cgl/users/mv2850/mixcr/mixcr

# Directory containing the input files
input_dir="/burg/cgl/users/mv2850/runs/250903_mv_jurkat_adept//mixcr"
output_dir="/burg/cgl/users/mv2850/runs/250903_mv_jurkat_adept/mixcr_barcodes/results"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Array to hold file prefixes
file_array=()
# Loop through the R1 files and store their prefixes in the array
for r1_file in "$input_dir"/*_R1.fastq.gz; do
    base=$(basename "$r1_file" _R1.fastq.gz)
    file_array+=("$base")
done

# Process in batches of 10
batch_size=10
for ((i=0; i<${#file_array[@]}; i+=batch_size)); do
    # Process each batch
    for ((j=i; j<i+batch_size && j<${#file_array[@]}; j++)); do
        base="${file_array[$j]}"
        echo "Processing: $base"  # Optional: print which base is being processed
        mixcr analyze local:sci-rna-seq3 \
            "$input_dir/${base}_R1.fastq.gz" \
            "$input_dir/${base}_R2.fastq.gz" \
            "$output_dir/TCR_sci_${base}"
    done
done

