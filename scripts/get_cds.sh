#!/bin/bash

# Check for user input
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input_directory> <api-key>"
    exit 1
fi

input_dir="$1"
api_key="$2"

# Validate input directory
if [[ ! -d "$input_dir" ]]; then
    echo "Error: '$input_dir' is not a valid directory."
    exit 1
fi

# Get all unique chunk IDs (basename without extension)
chunk_ids=$(find "$input_dir" -type f -name 'chunk_*' \
    | sed -E 's/\.zip$//' \
    | sed -E 's/\..*$//' \
    | xargs -n1 basename \
    | sort -u)

# Loop over each unique chunk ID
for chunk_id in $chunk_ids; do
    chunk_path="$input_dir/$chunk_id"
    zip_file="${chunk_path}.zip"

    if [[ -e "$zip_file" ]]; then
        if unzip -t "$zip_file" > /dev/null 2>&1; then
            echo "Skipping valid zip: $zip_file"
            continue
        else
            echo "Removing corrupted zip: $zip_file"
            rm -f "$zip_file"
        fi
    else
        echo "Downloading: $chunk_id"
    fi

    # Download if zip is missing or was removed
    datasets download genome accession \
        --inputfile "$chunk_path" \
        --include cds \
        --api-key "${api_key}" \
        --filename "$zip_file"

    echo "Finished downloading $chunk_id"
done
