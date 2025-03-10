#!/bin/bash
#==========================================================================
# RNA-Seq Analysis Pipeline for Streptococcus CHX Project
#==========================================================================
# Downloads and processes RNA-seq data (PRJNA1162077)
# Includes download, QC, trimming, annotation and differential expression
#
# Usage: ./download_and_process.sh [OPTIONS]
#
# Options:
#   --skip-download     Skip the download step
#   --skip-trimming     Skip the trimming step
#   --skip-gff          Skip the GFF preprocessing step
#   --strains ID1,ID2   Process only specified strains (comma-separated)
#   --help              Display this help message
#
#==========================================================================

# Set strict error handling
set -e
set -o pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Logging function
log() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

# Error handling function
error_exit() {
    log "ERROR: $1" >&2
    exit "${2:-1}"
}

# Help function
show_help() {
    grep '^#' "$0" | grep -v '^#!/bin/bash' | sed 's/^# //' | sed 's/^#//'
    exit 0
}

# Parse command line options
SKIP_DOWNLOAD=false
SKIP_TRIMMING=false
SKIP_GFF=false
STRAINS="73 78 93"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --skip-download) SKIP_DOWNLOAD=true ;;
        --skip-trimming) SKIP_TRIMMING=true ;;
        --skip-gff) SKIP_GFF=true ;;
        --strains) STRAINS=$(echo "$2" | tr ',' ' '); shift ;;
        --help) show_help ;;
        *) error_exit "Unknown parameter: $1" ;;
    esac
    shift
done

# Validate required directories
[ -d "$SCRIPT_DIR/Data" ] || mkdir -p "$SCRIPT_DIR/Data" || error_exit "Cannot create Data directory"

# Step 1: Download RNA-seq data
if [ "$SKIP_DOWNLOAD" = false ]; then
        log "Starting RNAseq data download..."
        bash "$SCRIPT_DIR/bash/download_rnaseq_data.sh" || error_exit "RNAseq download failed"
        log "RNAseq data download completed."
else
    log "Skipping download step as requested."
fi

# Step 2: Quality Control and Trimming
if [ "$SKIP_TRIMMING" = false ]; then
    if [ ! -d "Data/PRJNA1162077_data/trimmed" ]; then
        log "Trimming..."
        bash "$SCRIPT_DIR/bash/trim_adapters.sh" -i "Data/PRJNA1162077_data" -t $(( $(nproc) - 1 )) || error_exit "Cutadapt trim failed"
        log "Cutadapt trimming completed."
    else
        log "Trimmed files already exist, skipping cutadapt step."
    fi
else
    log "Skipping trimming step as requested."
fi

# Step 3: GFF File Preprocessing
if [ "$SKIP_GFF" = false ]; then
    for strain_id in $STRAINS; do
        output_dir="Data/Bakta_annotated_files/Strain_${strain_id}/Processed_for_nfcoreRnaseq"
        if [ ! -d "$output_dir" ]; then
            log "Preprocessing bacterial GFF files for Strain_${strain_id}..."
            bash "$SCRIPT_DIR/bash/preprocess_gff.sh" -i "Data/Bakta_annotated_files/Strain_${strain_id}" || 
                error_exit "Preprocessing GFF for Strain_${strain_id} failed"
        else
            log "Processed files already exist for Strain_${strain_id}, skipping preprocessing."
        fi
    done
    log "Bacterial GFF preprocessing completed for all selected strains."
else
    log "Skipping GFF preprocessing step as requested."
fi

# Step 4: Execute nf-core/rnaseq Pipeline
for strain_id in $STRAINS; do
    results_dir="results_strain_${strain_id}"
    if [ ! -d "$results_dir" ]; then
        log "Starting nf-core/rnaseq pipeline for Strain_${strain_id}..."
        bash "$SCRIPT_DIR/bash/run_nfcore_rnaseq.sh" "${strain_id}" --resume || 
            error_exit "Nextflow pipeline run failed for Strain_${strain_id}"
    else
        log "Results already exist for Strain_${strain_id}, skipping pipeline run."
    fi
done
log "nf-core/rnaseq pipeline completed successfully for all selected strains."

log "Preprocessing Successful. Execute the Main_data_analysis.R script for further analysis."

