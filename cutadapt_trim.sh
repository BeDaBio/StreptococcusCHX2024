#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FOLDER [-h]"
    echo "  -i INPUT_FOLDER  Specify the input folder containing all fastq.gz files"
    echo "  -h               Display this help message"
    exit 1
}

# Parse command-line options
while getopts ":i:h" opt; do
    case $opt in
        i) INPUT_FOLDER="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if INPUT_FOLDER is provided
if [ -z "$INPUT_FOLDER" ]; then
    echo "Error: Input folder not specified."
    usage
fi

# Check if cutadapt is available
if ! command -v cutadapt &> /dev/null; then
    echo "Error: cutadapt is not installed or not in the system PATH."
    exit 1
fi

# Define the output folder and create it if it doesn't exist
OUTPUT_FOLDER="$INPUT_FOLDER/trimmed"
mkdir -p "$OUTPUT_FOLDER"

# Loop through all input files and run Cutadapt
for input_file in "$INPUT_FOLDER"/*.fastq.gz; do
    # Check if there are any matching files
    if [ ! -e "$input_file" ]; then
        echo "Error: No .fastq.gz files found in $INPUT_FOLDER"
        exit 1
    fi

    # Extract the base name of the file (without the folder and extension)
    base_name=$(basename "$input_file" .fastq.gz)
  
    # Define the output file path
    output_file="$OUTPUT_FOLDER/${base_name}.AD2trimmed.fastq.gz"

    # Run Cutadapt with the specified parameters
    cutadapt -a CCATCTCATCCCTGCGTGTCTCCGACTCAG \
        -o "$output_file" \
        "$input_file" \
        -j 0

    echo "Trimmed file saved as $output_file"
done

echo "Cutadapt trimming completed for all files in $INPUT_FOLDER."


