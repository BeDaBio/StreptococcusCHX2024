# Main Script for Bioinformatics Analysis
# This script orchestrates our entire workflow - it loads organism data,
# creates Gene Ontology databases, and produces differential expression reports.
# The results are stored in the "output" folder.

library(here)
# Define relevant paths ---------------------------------------------------
# Path to the nf-core output directory
nf_core_output_path <- file.path(here(), "nf_core_output")
# Expected structure: nf_core_output/Strain_73, nf_core_output/Strain_78, nf_core_output/Strain_93

# Load required libraries -------------------------------------------------
# For consistent file paths across different operating systems
source(file.path(here(), "R/Libraries.R"))
source(file.path(here(), "R/CreateOrgDb.R"))

# Define variables --------------------------------------------------------
# Define organisms with their strain names and taxonomic IDs

# Named vector linking organism codes to strain names
organisms <- c(
  "73" = "Streptococcus salivariusS73wt",
  "78" = "Streptococcus vestibularisS78wt",
  "93" = "Streptococcus mitisS93wt"
)

# Extract strain numbers
strains <- gsub("[^0-9]", "", organisms)

# Named vector linking organism codes to NCBI Taxonomy IDs
tax_ids <- c("73" = 1304, "78" = 1343, "93" = 28037)

# Set up directories -----------------------------------------------------
# Create directory structure for the project

# GO annotation directories
obo_dir <- file.path(here(), "Data", "GO_Check")
go_obo_file <- file.path(obo_dir, "go-basic.obo")
do_not_annotate_file <- file.path(obo_dir, "gocheck_do_not_annotate.obo")

# Create directories if they don't exist
if (!dir.exists(obo_dir)) {
  dir.create(obo_dir, recursive = TRUE)
}

# Output directories
output_quarto <- file.path(here(), "output", "Quarto_Reports")
if (!dir.exists(output_quarto)) {
  dir.create(output_quarto, recursive = TRUE)
}

# Organism databases path
org_db_path <- file.path(here(), "OrgDbs")
if (!dir.exists(org_db_path)) {
  dir.create(org_db_path, recursive = TRUE)
}

# Download GO ontology files ---------------------------------------------
# Download files if they don't exist

# Download gocheck_do_not_annotate.obo
if (!file.exists(do_not_annotate_file)) {
  download.file(
    "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.obo",
    destfile = do_not_annotate_file
  )
}

# Download go-basic.obo
if (!file.exists(go_obo_file)) {
  download.file(
    "https://purl.obolibrary.org/obo/go/go-basic.obo",
    destfile = go_obo_file
  )
}

# Load or create organism databases --------------------------------------
# Check if organism databases exist and load them, otherwise create new ones

# List existing organism databases
org_db_files <- list.files(org_db_path, pattern = "org.*.db", full.names = FALSE)

# Required database names
required_dbs <- c(
  "org.SmitisS93wt.eg.db",
  "org.SsalivariusS73wt.eg.db",
  "org.SvestibularisS78wt.eg.db"
)

# If all databases exist, load them
if (all(required_dbs %in% org_db_files)) {
  for (db in org_db_files) {
    library(db, character.only = TRUE)
  }
} else {
  # Create missing databases
  mapply(function(tax_id, organism) {
    # Generate strain ID
    strain <- paste(gsub("[^0-9]", "", organism), "wt", sep = "_")

    # Find GO annotation file
    pannzer_files <- list.files(
      file.path(here(), "Data", "GO_Pannzer2"),
      pattern = strain,
      recursive = TRUE,
      full.names = TRUE
    )
    pannzer_file <- grep("GO_strain", pannzer_files, value = TRUE)

    # If no GO file found, set to NA
    if (length(pannzer_file) == 0) {
      pannzer_file <- NA
    }

    # Read Bakta annotation file
    bakta_file <- list.files(
      file.path(here(), "Data", "Bakta_annotated_files"),
      pattern = paste0(strain, ".tsv"),
      recursive = TRUE,
      full.names = TRUE
    )

    bakta_table <- read_tsv(bakta_file, skip = 5)

    # Clean column names
    colnames(bakta_table) <- gsub("#", "", colnames(bakta_table))
    colnames(bakta_table) <- gsub(" ", "_", colnames(bakta_table))

    # Create organism database
    createOrgfromBakta(
      outdir = org_db_path,
      Strain = strain,
      Pannzer = pannzer_file,
      Organism = organism,
      Bakta_tsv = bakta_table,
      GID = "Locus_Tag",
      EVIDENCE = "IEA",
      GOID = "GO",
      install_DB = FALSE,
      tax_id = tax_id
    )
  }, tax_ids, organisms)

  # Install newly created databases
  org_db_paths <- list.files(file.path(here(), "OrgDbs"), pattern = "org.*.db", full.names = TRUE)
  install.packages(org_db_paths, repos = NULL)

  # Load the databases
  lapply(basename(org_db_paths), require, character.only = TRUE)
}

# Differential expression analysis ---------------------------------------
# Generate differential expression reports using Quarto

process_quarto_render <- function(input, output_dir = NULL, output_file = NULL, ...) {
  # Inspect the Quarto file
  x <- quarto::quarto_inspect(input)
   output_format <- names(x$formats)
  # # Determine output file name
  # output_file <- if (is.null(output_file)) {
  #   
  # } else {
  #   paste0(output_file,".",output_format)
  # }
  # print(output_file)
  # Render the document
  quarto::quarto_render(input = input, ... )#, output_file = output_file,execute_dir = file.path(here(),"R")
  # Determine input and output directories
  input_dir <- dirname(input)
  output_dir <- if (is.null(output_dir)) input_dir else output_dir
  
  # Handle file location
  output_from <- file.path("R",x$formats[[output_format]]$pandoc$`output-file`)
  output_to <- file.path(output_dir, paste0(output_file,".",output_format))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.copy(output_from, output_to, overwrite = TRUE)) {
    file.remove(output_from)
  }

}

# Render reports for each strain
for (strain in strains) {
  process_quarto_render(
    input = "R/DifferentialExpression.qmd",
    output_file = paste0("DifferentialExpression_Strain", strain),
    output_dir = output_quarto,
    execute_params = list(
      "strain" = strain,
      "nfcore_output_Path" = nf_core_output_path
    )#,
   # execute_dir = file.path(here(),"R")
  )
}

# Cross-strain comparison ------------------------------------------------
# Generate cross-strain comparison report

process_quarto_render(
  input = "R/CrossStrain_Analyis.qmd",
  output_file = "Cross-strain_Analysis",
  output_dir = output_quarto,
  execute_params = list()
)
