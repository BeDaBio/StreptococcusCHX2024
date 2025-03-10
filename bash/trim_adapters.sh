#!/bin/bash
#===============================================================================
# RNA-Seq Adapters Trimming Script
#===============================================================================
# This script processes RNA-seq data files by removing adapter sequences and 
# performing quality trimming using Trimmomatic and Cutadapt.
#
# Usage: 
#   ./trim_adapters.sh [OPTIONS]
#
# Options:
#   -i, --input-dir DIR     Directory containing the fastq.gz files (required)
#   -o, --output-dir DIR    Directory to save trimmed data (default: input-dir/trimmed)
#   -t, --threads NUM       Number of threads to use (default: 4)
#   -h, --help              Display this help message
#
# Dependencies:
# - Trimmomatic
# - Cutadapt
#
#===============================================================================

set -e
set -o pipefail

# Default parameters
INPUT_FOLDER=""
OUTPUT_FOLDER=""
THREADS=4

# Function to display help message
show_help() {
    grep '^#' "$0" | grep -v '^#!/bin/bash' | sed 's/^# //' | sed 's/^#//'
    exit 0
}

# Function for logging
log() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

# Function for error handling
error_exit() {
    log "ERROR: $1" >&2
    exit "${2:-1}"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-dir)
            INPUT_FOLDER="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_FOLDER="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            error_exit "Unknown option: $1"
            ;;
    esac
done

# Check required parameters
if [ -z "$INPUT_FOLDER" ]; then
    error_exit "Input folder not specified. Use -i or --input-dir option."
fi

# Set default output folder if not specified
if [ -z "$OUTPUT_FOLDER" ]; then
    OUTPUT_FOLDER="$INPUT_FOLDER/trimmed"
fi

# Create output directory
mkdir -p "$OUTPUT_FOLDER" || error_exit "Cannot create output directory: $OUTPUT_FOLDER"

# Check required tools
for tool in cutadapt trimmomatic; do
    command -v $tool >/dev/null 2>&1 || error_exit "$tool not found. Please install required dependencies."
done

# Loop through all input files and run processing
log "Starting adapter trimming workflow..."
FILE_COUNT=0

for input_file in "$INPUT_FOLDER"/*.fastq.gz; do
    # Check if there are any matching files
    if [ ! -e "$input_file" ]; then
        error_exit "No .fastq.gz files found in $INPUT_FOLDER"
    fi

    # Extract the base name of the file
    SAMPLE=$(basename "$input_file" .fastq.gz)
    log "Processing $SAMPLE"
    
    # Uncompress input file for processing
    gunzip -c "$input_file" > "$OUTPUT_FOLDER/${SAMPLE}.fastq"

    # Step 1: Quality Trimming with Trimmomatic
    log "Running Trimmomatic quality trimming for $SAMPLE..."
    trimmomatic SE "$OUTPUT_FOLDER/${SAMPLE}.fastq" \
        "$OUTPUT_FOLDER/${SAMPLE}.trimmed.fastq" \
        -phred33 -threads $THREADS \
        SLIDINGWINDOW:25:15 MINLEN:80 \
        -trimlog "$OUTPUT_FOLDER/trimlog_${SAMPLE}.log"

    # Step 2: P1 Adapter Removal
    log "Running first adapter removal (P1) for $SAMPLE..."
    cutadapt -a CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT \
        -o "$OUTPUT_FOLDER/${SAMPLE}.AD1trimmed.fastq" \
        "$OUTPUT_FOLDER/${SAMPLE}.trimmed.fastq" \
        -j $THREADS > "$OUTPUT_FOLDER/Adapter_trimming_P1_${SAMPLE}.log"

    # Step 3: A Adapter Removal
    log "Running second adapter removal (A) for $SAMPLE..."
    cutadapt -a CCATCTCATCCCTGCGTGTCTCCGACTCAG \
        -o "$OUTPUT_FOLDER/${SAMPLE}.AD2trimmed.fastq" \
        "$OUTPUT_FOLDER/${SAMPLE}.AD1trimmed.fastq" \
        -j $THREADS > "$OUTPUT_FOLDER/Adapter_trimming_A_${SAMPLE}.log"

    # Compress final output
    log "Compressing final output for $SAMPLE..."
    gzip -f "$OUTPUT_FOLDER/${SAMPLE}.AD2trimmed.fastq"
    
    # Cleanup intermediate files
    rm -f "$OUTPUT_FOLDER/${SAMPLE}.fastq"
    rm -f "$OUTPUT_FOLDER/${SAMPLE}.trimmed.fastq"
    rm -f "$OUTPUT_FOLDER/${SAMPLE}.AD1trimmed.fastq"
    
    log "Completed processing $SAMPLE"
    FILE_COUNT=$((FILE_COUNT + 1))
done

log "Adapter trimming completed for $FILE_COUNT files."

# Update sample sheets to point to trimmed files
log "Updating sample sheet file paths..."
for sample_csv in "$INPUT_FOLDER"/samplefile_*.csv; do
    if [ -f "$sample_csv" ]; then
        # Update the fastq file paths in the CSV file
        # Extract the directory name without the full INPUT_FOLDER path
        RELATIVE_DIR=$(echo "$(dirname "$sample_csv")" | sed "s|^$INPUT_FOLDER/||")
        
        # Replace original fastq paths with trimmed ones
        sed -i -E "s#(Data/PRJNA1162077_data)/([^/]+\.fastq\.gz)#\1/trimmed/\2#g" "$sample_csv"
        sed -i -E "s#\.fastq\.gz#\.AD2trimmed\.fastq\.gz#g" "$sample_csv"
        log "Updated paths in $sample_csv"
    fi
done

log "All processing complete. Trimmed files are in $OUTPUT_FOLDER"