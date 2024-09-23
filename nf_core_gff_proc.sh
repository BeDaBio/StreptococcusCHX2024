#!/bin/bash

# Processing GTF/GFF/GFF3 files for Nextflow nf-core/rnaseq 3.14
# Display help
display_help() {
    echo "Usage: $0 -i <path to directory containing .gff/.gff3/.gtf and corresponding genome .fna file> -o <output directory>"
    echo "Make sure 'gffread' is installed and accessible in the system path."
    exit 1
}

# Parse command line arguments
while getopts "hi:o:" opt; do
  case ${opt} in
    h ) display_help ;;
    i ) input_dir=$OPTARG ;;
    o ) output_dir=$OPTARG ;;
    \? ) echo "Invalid option: -$OPTARG" 1>&2; display_help ;;
    : ) echo "Option -$OPTARG requires an argument." 1>&2; display_help ;;
  esac
done
shift $((OPTIND -1))

# Check if input_dir is provided
if [ -z "$input_dir" ]; then
    echo "Error: Input directory is required."
    display_help
fi

# Set default output_dir if not provided
output_dir=${output_dir:-"$input_dir/Processed_for_nfcoreRnaseq"}
mkdir -p "$output_dir"

# Ensure gffread is installed
if ! command -v gffread >/dev/null 2>&1; then
    echo "Error: gffread is not installed or not in your PATH."
    exit 1
fi

# Find input files
gtf_file=$(find "$input_dir" -type f -name "*.gtf" | head -n 1)
gff_file=$(find "$input_dir" -type f \( -name "*.gff" -o -name "*.gff3" \) | head -n 1)
genome_file=$(find "$input_dir" -type f -name "*.fna" | head -n 1)

if [ -z "$genome_file" ]; then
    echo "Genome file (*.fna) not found in the input directory."
    exit 1
fi

# Process GTF/GFF file
if [ -n "$gtf_file" ]; then
    input_file=$gtf_file
elif [ -n "$gff_file" ]; then
    input_file=$gff_file
    # Convert GFF/GFF3 to GTF using gffread (output to stdout)
    gtf_content=$(gffread -E -F -T "$input_file")
else
    echo "No GTF, GFF, or GFF3 file found in the input directory."
    exit 1
fi

# Filter for entries with transcript_id and save as the final GTF
final_gtf="$output_dir/filtered_transcripts.gtf"
if [ -n "$gtf_file" ]; then
    awk '$9 ~ /transcript_id/ {print}' "$input_file" > "$final_gtf"
else
    echo "$gtf_content" | awk '$9 ~ /transcript_id/ {print}' > "$final_gtf"
fi

# Generate a transcript FASTA file using gffread
transcript_fasta="$output_dir/transcripts.fna"
gffread -w "$transcript_fasta" -g "$genome_file" "$final_gtf"

echo "Processing completed. Output files are in $output_dir"
