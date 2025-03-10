#!/bin/bash
#===============================================================================
# RNA-Seq Data Download Script
#===============================================================================
# This script downloads RNA-seq data from NCBI SRA database for project PRJNA1162077
# and organizes the files in the appropriate directory structure.
#
# Usage: 
#   ./download_RNAseq_data.sh [OPTIONS]
#
# Options:
#   -o, --output-dir DIR    Directory to save downloaded data (default: Data/PRJNA1162077_data)
#   -t, --threads NUM       Number of threads to use for download (default: 1)
#   -h, --help              Display this help message
#
# Dependencies:
# - SRA Toolkit (prefetch, fasterq-dump)
# - curl or wget
#
# Author: Your Name
# Date: YYYY-MM-DD
#===============================================================================

set -e
set -o pipefail

# Default parameters
OUTPUT_DIR="Data/PRJNA1162077_data"
THREADS=1
ACCESSION="PRJNA1162077"
echo "$OUTPUT_DIR/${SRR_ID}.fastq.gz"
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
        -o|--output-dir)
            OUTPUT_DIR="$2"
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

# Validate SRA tools existence
command -v prefetch >/dev/null 2>&1 || error_exit "SRA Toolkit's prefetch command not found. Please install SRA Toolkit."
command -v fasterq-dump >/dev/null 2>&1 || error_exit "SRA Toolkit's fasterq-dump command not found. Please install SRA Toolkit."

# Create output directory
mkdir -p "$OUTPUT_DIR" || error_exit "Cannot create output directory: $OUTPUT_DIR"

# Download SRA run information
log "Fetching run information for $ACCESSION..."
RUN_INFO_FILE="$OUTPUT_DIR/runinfo.csv"
esearch -db sra -query "$ACCESSION" | efetch -format runinfo > "$RUN_INFO_FILE" || 
    error_exit "Failed to fetch run information for $ACCESSION"

# Extract SRR IDs
SRR_IDS=$(tail -n +2 "$RUN_INFO_FILE" | cut -d',' -f1)
SRR_COUNT=$(echo "$SRR_IDS" | wc -l)
log "Found $SRR_COUNT runs to download"

# Download and convert SRA files to FASTQ
for SRR_ID in $SRR_IDS; do
    if [[ ! -f "$OUTPUT_DIR/${SRR_ID}.fastq.gz" ]]; then
        log "Downloading $SRR_ID..."
        
        # Create a temporary directory for the download
        TMP_DIR=$(mktemp -d)
        
        # Download SRA file
        prefetch "$SRR_ID" -O "$TMP_DIR" || error_exit "Failed to download $SRR_ID"
        
        # Convert to FASTQ
        log "Converting $SRR_ID to FASTQ..."
        fasterq-dump --split-files --threads "$THREADS" --outdir "$TMP_DIR" "$TMP_DIR/$SRR_ID/$SRR_ID.sra" || 
            error_exit "Failed to convert $SRR_ID to FASTQ"
        
        # Compress FASTQ files
        log "Compressing FASTQ files for $SRR_ID..."
        for FASTQ_FILE in "$TMP_DIR/"*.fastq; do
            gzip -c "$FASTQ_FILE" > "$OUTPUT_DIR/$(basename "$FASTQ_FILE").gz" || 
                error_exit "Failed to compress $FASTQ_FILE"
        done
        
        # Clean up temporary files
        rm -rf "$TMP_DIR"
        
        log "Successfully processed $SRR_ID"
    else
        log "$SRR_ID already downloaded, skipping"
    fi
done

# Create metadata file with experiment information
log "Creating metadata file..."
# The BioSample column is at the 26th place in runinfo.csv
awk -F ',' 'NR>1 {print $26}' "$OUTPUT_DIR/runinfo.csv" > "$OUTPUT_DIR/biosample_ids.txt"

# Write header with the renamed Treatment column (8 columns) as CSV
echo "BioSample,SRR,Treatment,Organism,Strain,Genotype,CHX,Rep" > "$OUTPUT_DIR/biosample_metadata.csv"

while read -r BS; do
  # Skip empty lines
  if [[ -n "$BS" ]]; then
    echo "Processing $BS"
    # Create a temporary file to preserve multi-line output
    TMPFILE="$(mktemp)"

    SRR="$(awk -F ',' -v bs="$BS" 'NR>1 && $26==bs {print $1; exit}' "$OUTPUT_DIR/runinfo.csv")"

    esearch -db biosample -query "$BS" < /dev/null 2>/dev/null \
      | efetch -format docsum 2>/dev/null \
      | xtract \
         -pattern DocumentSummary \
         -element Title Organism \
         -block Attribute \
           -if Attribute@attribute_name -equals "strain" \
           -element Attribute \
      > "$TMPFILE"

    # Check if we got any output
    if [[ -s "$TMPFILE" ]]; then
      # Get Treatment value from the first column of TMPFILE
      Treatment=$(head -1 "$TMPFILE" | awk '{print $1}')
      echo $Treatment
      # Retrieve organism information
      Organism=$(grep -o "Streptococcus.*" "$TMPFILE" | awk '{print $1, $2}' || echo "Unknown")
      echo $Organism
      # Split Treatment to extract its components
      IFS='_' read -r -a parts <<< "$Treatment"
      Strain="${parts[0]:-}"
      Genotype="${parts[1]:-}"
      CHX="${parts[2]:-}"
      Rep="${parts[3]:-}"
      
      # Output metadata with a single Strain column and renamed Treatment column
      echo -e "$BS,$SRR,$Treatment,$Organism,$Strain,$Genotype,$CHX,$Rep" >> "$OUTPUT_DIR/biosample_metadata.csv"
    else
      # Use commas in the "No record found" line
      echo -e "$BS,$SRR,No record found,,,,," >> "$OUTPUT_DIR/biosample_metadata.csv"
    fi

    rm -f "$TMPFILE"
  fi
done < "$OUTPUT_DIR/biosample_ids.txt"

echo "Creating strain-specific sample files..."
# Use awk to extract the strain column (column 5) from CSV
awk -F',' 'NR>1 {print $5}' "$OUTPUT_DIR/biosample_metadata.csv" | sort -u | while read -r strain; do
  if [[ -n "$strain" ]]; then
    outfile="$OUTPUT_DIR/samplefile_${strain}.csv"
    # Write header for sample files
    echo "sample,fastq_1,fastq_2,strandedness,BioSample,Strain,Genotype,CHX,Rep,Organism" > "$outfile"
    # Extract matching lines and parse Treatment into components
    awk -v outdir="$OUTPUT_DIR" -v strain="$strain" -F',' '
      $5 == strain {
        # Mapping: $1 BioSample, $2 SRR, $3 Treatment, $4 Organism,
        # $5 Strain, $6 Genotype, $7 CHX, $8 Rep
        printf "%s,%s/%s.fastq.gz,,%s,%s,%s,%s,%s,%s\n", $2, outdir, $2, "auto", $1, $5, $6, $7, $8, $4
      }' "$OUTPUT_DIR/biosample_metadata.csv" >> "$outfile"
    echo "Created $outfile"
  fi
done


log "Download complete. Data saved to $OUTPUT_DIR"
log "Total files downloaded: $SRR_COUNT"

# Validate the download
EXPECTED_FILES=$((SRR_COUNT)) 
ACTUAL_FILES=$(find "$OUTPUT_DIR" -name "*.fastq.gz" | wc -l)
if [[ "$ACTUAL_FILES" -eq "$EXPECTED_FILES" ]]; then
    log "Download validation successful: All expected files were downloaded"
else
    log "WARNING: Expected $EXPECTED_FILES files but found $ACTUAL_FILES. Some files may be missing."
fi
