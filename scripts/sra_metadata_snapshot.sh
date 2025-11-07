#! /bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_csv_file>"
    exit 1
fi

echo 
input_csv="$1"

if [ ! -f "$input_csv" ]; then
    echo "Error: File '$input_csv' not found"
    exit 1
fi

# SRA metadata can have additional commas included between '()'
# replace the comma with an underscore
header=$(awk 'NR==1 {
    while (match($0, /\([^()]*,[^()]*\)/)) {
        s = substr($0, RSTART, RLENGTH);
        gsub(",", "_", s);
        $0 = substr($0, 1, RSTART - 1) s substr($0, RSTART + RLENGTH);
    }
    print;
    next
} 1' $input_csv)

IFS=',' read -ra columns <<< "$header"

num_cols=${#columns[@]}
echo "Found '$num_cols' columns in '$input_csv' file"

#https://stackoverflow.com/a/31889595
FPAT='([^,]*)|("[^"]*")'

for ((i=1; i<=num_cols; i++)); do
    col_name="${columns[i-1]}"
#    if [[ "$col_name" == "Run" ]]; then
#        continue
#    fi
    unique_count=$(awk -v FPAT="$FPAT" -v col="$i" 'NR > 1 {if ($col != "") print $col}' "$input_csv" | sort | uniq | wc -l)
    total_count=$(awk -v FPAT="$FPAT" -v col="$i" 'NR > 1 {if ($col != "") print $col}' "$input_csv" | wc -l)
    echo "=== $unique_count Unique values (of $total_count values) for column $i: $col_name ==="
#    awk -F',' -v col="$i" 'NR > 1 {print $col}' "$input_csv" | sort | uniq -c | sort -nr
    awk -v FPAT="$FPAT" -v col="$i" 'NR > 1 {if ($col != "") print $col}' "$input_csv" | sort | uniq -c | sort -nr
    echo
done
