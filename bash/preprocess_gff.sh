#!/bin/bash
#===============================================================================
# Bacterial GFF Preprocessing Script
#===============================================================================
# Preprocesses bacterial GFF files for compatibility with nf-core/rnaseq pipeline
# Handles:
# - Gene/transcript ID matching
# - Missing transcript features
# - Feature hierarchy
# - Biotype conversion
#
# Usage: 
#   ./preprocess-gff.sh [OPTIONS]
#
# Options:
#   -i, --input-dir DIR       Directory with input Bakta annotation files (required)
#   -o, --output-dir DIR      Output directory (default: input-dir/Processed_for_nfcoreRnaseq)
#   -g, --genome-fasta FILE   Path to genome FASTA file (default: auto-detect)
#   -h, --help                Display this help message
#
# Dependencies: gffread, bedtools, awk, sed
#===============================================================================

set -e
set -o pipefail

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
GENOME_FASTA=""
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
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--genome-fasta)
            GENOME_FASTA="$2"
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


echo $(pwd)
# Validate input parameters
[[ -z "$INPUT_DIR" ]] && error_exit "Input directory is required (-i, --input-dir)"
[[ ! -d "$INPUT_DIR" ]] && error_exit "Input directory does not exist: $INPUT_DIR"

# Set default output directory if not specified
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${INPUT_DIR}/Processed_for_nfcoreRnaseq"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR" || error_exit "Cannot create output directory: $OUTPUT_DIR"

# Validate dependency installations
command -v gffread >/dev/null 2>&1 || error_exit "gffread not found. Please install gffread."
command -v bedtools >/dev/null 2>&1 || error_exit "bedtools not found. Please install bedtools."

# Find input GFF file
GFF_FILE=$(find "$INPUT_DIR" -name "*.gff3" -o -name "*.gff" | head -n 1)
[[ -z "$GFF_FILE" ]] && error_exit "No GFF file found in input directory: $INPUT_DIR"
log "Found GFF file: $GFF_FILE"

# Autodetect genome FASTA file if not provided
if [[ -z "$GENOME_FASTA" ]]; then
    GENOME_FASTA=$(find "$INPUT_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" | head -n 1)
    [[ -z "$GENOME_FASTA" ]] && error_exit "No genome FASTA file found. Please specify with -g, --genome-fasta"
    log "Found genome FASTA file: $GENOME_FASTA"
fi

# Output files
PROCESSED_GFF="${OUTPUT_DIR}/processed_annotation.gtf"
GENOME_OUTPUT="${OUTPUT_DIR}/genome.fasta"

# Step 1: Copy and standardize genome FASTA file
log "Copying and standardizing genome FASTA file..."
cat "$GENOME_FASTA" | awk '/^>/ {gsub(/[^a-zA-Z0-9_.]/, "_", $0); print $0} /^[^>]/ {print}' > "$GENOME_OUTPUT"

# Step 2: Convert GFF3 to GTF with necessary modifications
log "Converting GFF3 to GTF format with modifications..."
TMP_GFF="${OUTPUT_DIR}/tmp.gff"

# First, add transcript features if they are missing
awk '
BEGIN {OFS="\t"}
{
    if ($3 == "gene" || $3 == "CDS" || $3 == "exon") {
        print $0
        if ($3 == "gene") {
            # Extract gene ID
            match($9, /ID=([^;]+)/)
            if (RSTART > 0) {
                gene_id = substr($9, RSTART+3, RLENGTH-3)
                # Create a transcript feature with same coordinates
                transcript_id = gene_id ".t"
                attributes = "ID=" transcript_id ";Parent=" gene_id
                print $1, $2, "transcript", $4, $5, $6, $7, $8, attributes
            }
        }
    } else {
        print $0
    }
}' "$GFF_FILE" > "$TMP_GFF"

# Convert to GTF using gffread
gffread -E "$TMP_GFF" -T -o "$PROCESSED_GFF" || error_exit "Failed to convert GFF to GTF"

# Further process the GTF file to ensure compatibility
TMP_GTF="${OUTPUT_DIR}/tmp.gtf"
awk '
{
    if ($3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS") {
        # Find gene_id
        match($0, /gene_id "([^"]+)"/)
        if (RSTART > 0) {
            gene_id = substr($0, RSTART+9, RLENGTH-10)
            
            # Find or create transcript_id
            match($0, /transcript_id "([^"]+)"/)
            if (RSTART > 0) {
                transcript_id = substr($0, RSTART+15, RLENGTH-16)
            } else {
                transcript_id = gene_id ".t"
            }
            
            # Add biotype for gene and transcript
            biotype = "protein_coding"
            
            # Rebuild attributes
            attrs = " gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\"; "
            
            if ($3 == "gene") {
                attrs = attrs "gene_biotype \"" biotype "\"; "
            } else if ($3 == "transcript") {
                attrs = attrs "transcript_biotype \"" biotype "\"; "
            }
            
            # Print modified line
            print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" attrs
        } else {
            print
        }
    } else {
        print
    }
}' "$PROCESSED_GFF" > "$TMP_GTF"

mv "$TMP_GTF" "$PROCESSED_GFF"

# Clean up temporary files
rm -f "$TMP_GFF"

# Create gene-level table
log "Creating gene information table..."
GENE_TABLE="${OUTPUT_DIR}/gene_info.tsv"

awk '$3 == "gene" {
    # Extract gene_id
    match($0, /gene_id "([^"]+)"/)
    if (RSTART > 0) {
        gene_id = substr($0, RSTART+9, RLENGTH-10)
        
        # Extract gene name if available
        gene_name = gene_id
        match($0, /gene_name "([^"]+)"/)
        if (RSTART > 0) {
            gene_name = substr($0, RSTART+11, RLENGTH-12)
        }
        
        # Extract biotype
        biotype = "protein_coding"
        match($0, /gene_biotype "([^"]+)"/)
        if (RSTART > 0) {
            biotype = substr($0, RSTART+14, RLENGTH-15)
        }
        
        print gene_id "\t" gene_name "\t" biotype
    }
}' "$PROCESSED_GFF" > "$GENE_TABLE"

log "Processing complete. Output files:"
log " - Processed GTF: $PROCESSED_GFF"
log " - Genome FASTA: $GENOME_OUTPUT"
log " - Gene Info Table: $GENE_TABLE"
