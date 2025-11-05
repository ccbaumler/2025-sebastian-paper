#!/bin/bash

input_csv="../metadata/SraRunTable-crc.txt"

# Define the patterns and corresponding file-safe names
declare -A categories=(
    ["healthy control"]="healthy_control.csv"
    ["colorectal cancer"]="colorectal_cancer.csv"
    ["Normal"]="Normal.csv"
    ["Cancer"]="Cancer.csv"
)

# Create each file with header
for category in "${!categories[@]}"; do
    cleanCategory=$(echo "$category"  | sed 's/ /_/g')
    echo -e "$cleanCategory" > "${categories[$category]}"
done

FPAT='([^,]*)|("[^"]*")'

for category in "${!categories[@]}"; do
    awk -v FPAT="$FPAT" -v col1=1 -v col57=57 -v target="$category" \
        'NR > 1 {
            gsub(/^ *| *$/, "", $col57)
            if ($col57 == target && $col1 != "")
                print $col1
        }' "$input_csv" >> "${categories[$category]}"
done
