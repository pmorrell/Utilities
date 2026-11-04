#!/bin/bash

# Path to the file containing the list of filenames
file_list="filenames.txt"

# Your S3 prefix string
s3_prefix="s3://morrelllab/NAM_Parents/"


# Loop through each line safely
while IFS= read -r file_name; do
    # Skip empty lines
    [[ -z "$file_name" ]] && continue
    
    echo "Processing file: $file_name"
  
      # Combine the prefix and the filename
    s3_path="${s3_prefix}${file_name}"
    
    s3cmd rm "$s3_path"

    
    # Your file processing logic goes here
    # Example: if [ -f "$file_name" ]; then ... fi

done < "$file_list"

