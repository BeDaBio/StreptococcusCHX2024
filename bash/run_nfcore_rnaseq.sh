#!/bin/bash
#===============================================================================
# nf-core/rnaseq Pipeline Runner
#===============================================================================
# This script executes the nf-core/rnaseq pipeline for bacterial RNA-seq data.
# It automatically configures the pipeline for a specific strain by setting
# appropriate parameters and input files.
#
# Usage: 
#   ./Nextflow.sh STRAIN_ID [OPTIONS]
#
# Arguments:
#   STRAIN_ID                Strain identifier (e.g., 73, 78, 93) (required)
#
# Options:
#   --memory MEM             Memory limit (default: 32.GB)
#   --cpus CPUS              Number of CPUs (default: 8)
#   --outdir DIR             Output directory (default: nf_core_output/Strain_STRAIN_ID)
#   --resume                 Resume previous run (default: false)
#   --help                   Display this help message
#
# Dependencies:
# - Nextflow
# - nf-core/rnaseq pipeline (version 3.14.0)
# - Singularity
#
#===============================================================================

set -e
set -o pipefail

# Default parameters
MEMORY="32.GB"
CPUS=8
RESUME=""
CUSTOM_OUTDIR=""

# Fixed parameters
PROFILE="singularity"
NFCORE_VERSION="3.14.0"

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
if [ $# -eq 0 ]; then
    show_help
fi

STRAIN="$1"
shift

# Validate strain parameter
if ! [[ "$STRAIN" =~ ^[0-9]+$ ]]; then
    error_exit "Invalid strain ID: $STRAIN. Must be a number."
fi

# Parse additional arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --memory)
            MEMORY="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --outdir)
            CUSTOM_OUTDIR="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        --help)
            show_help
            ;;
        *)
            error_exit "Unknown option: $1"
            ;;
    esac
done

# Set output directory
OUTDIR=${CUSTOM_OUTDIR:-"nf_core_output/Strain_${STRAIN}"}

# Validate required files
INPUT_FILE="Data/PRJNA1162077_data/samplefile_${STRAIN}.csv"
FASTA_FILE="Data/Bakta_annotated_files/Strain_${STRAIN}/Strain_${STRAIN}_wt.fna"
GTF_FILE="Data/Bakta_annotated_files/Strain_${STRAIN}/Processed_for_nfcoreRnaseq/Strain_${STRAIN}_wt_filtered.gtf"
TRANSCRIPT_FASTA="Data/Bakta_annotated_files/Strain_${STRAIN}/Processed_for_nfcoreRnaseq/Strain_${STRAIN}_wt_transcript.fna"

[[ ! -f "$INPUT_FILE" ]] && error_exit "Sample file not found: $INPUT_FILE"
[[ ! -f "$FASTA_FILE" ]] && error_exit "Genome FASTA file not found: $FASTA_FILE"
[[ ! -f "$GTF_FILE" ]] && error_exit "GTF file not found: $GTF_FILE"
[[ ! -f "$TRANSCRIPT_FASTA" ]] && error_exit "Transcript FASTA file not found: $TRANSCRIPT_FASTA"

# Validate nextflow installation
command -v nextflow >/dev/null 2>&1 || error_exit "Nextflow not found. Please install Nextflow."

# Create work directory for nextflow if it doesn't exist
WORK_DIR="work_strain_${STRAIN}"
mkdir -p "$WORK_DIR"

# Log the pipeline parameters
log "Running nf-core/rnaseq for strain $STRAIN"
log "Parameters:"
log " - Output directory: $OUTDIR"
log " - Memory: $MEMORY"
log " - CPUs: $CPUS" 
log " - Using singularity profile and nf-core/rnaseq version 3.14.0"
log " - Resume: ${RESUME:+Yes}"

# Run the nextflow command
log "Starting nextflow pipeline..."
nextflow run nf-core/rnaseq -r "$NFCORE_VERSION" \
  -profile "$PROFILE" \
  -work-dir "$WORK_DIR" \
  --input "$INPUT_FILE" \
  --outdir "$OUTDIR" \
  --fasta "$FASTA_FILE" \
  --gtf "$GTF_FILE" \
  --skip_umi_extract \
  --skip_trimming \
  --aligner star_salmon \
  --extra_star_align_args '--sjdbGTFfeatureExon CDS' \
  --featurecounts_feature_type transcript \
  --transcript_fasta "$TRANSCRIPT_FASTA" \
  --skip_dupradar \
  --skip_qualimap \
  --skip_rseqc \
  --featurecounts_group_type gene_id \
  --max_memory "$MEMORY" \
  --max_cpus "$CPUS" \
  $RESUME || error_exit "Nextflow pipeline failed"

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    log "Pipeline completed successfully for strain $STRAIN"
else
    error_exit "Pipeline failed for strain $STRAIN"
fi