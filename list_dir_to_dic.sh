#!/bin/bash

# Parse command-line arguments
while getopts ":d:o:" opt; do
  case $opt in
    d)
      directory=$OPTARG
      ;;
    o)
      output_file=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if directory and output file are provided
if [[ -z $directory || -z $output_file ]]; then
  echo "Usage: bash script.sh -d directory_to_list -o output.txt"
  exit 1
fi

# Change directory
cd "$directory" || exit

# Initialize an associative array
declare -A file_dict

# Iterate over each file in the directory
for file in *; do
    # Extract pattern from filename
    pattern=$(echo "$file" | awk -F '_' '{print $1"_"$2"_"$3}')
    # Append full path of the file to array element corresponding to its pattern
    file_dict["$pattern"]+="$directory$file "
done

# Write the pattern and associated files (with full paths) to the output file
for pattern in "${!file_dict[@]}"; do
    echo "$pattern,${file_dict[$pattern]}"
done > "$output_file"

echo "Files listed in pattern and saved to $output_file."
