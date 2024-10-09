# Main Script for Bioinformatics Analysis
# This script runs the entire bioinformatics pipeline, including loading organisms' data, 
# generating Gene Ontology databases, and producing differential expression reports.
# The results are stored in the "output" folder.

# Define relevant Paths ---------------------------------------------------

# nfcore_output_Path is the path to the nf-core output directory
nfcore_output_Path <- "nfcore_output" # Change this to the actual path where nf-core output is stored
# Best in a subfolder of the project directory with different subfolders for each strain 
# e.g. nfcore_output/Strain_73, nfcore_output/Strain_78, nfcore_output/Strain_93
# The Strain number in the subdirectory is used to grep the right data for each organism 

# Load Required Libraries --------------------------------------------------
library(here)  # Here is used for consistent file paths across different operating systems
source(file.path(here(), "Libraries.R"))  # Source a custom script with additional library dependencies or functions
source(file.path(here(), "CreateOrgDb.R"))  # Source a custom script with additional library dependencies or functions

# Define Variables ---------------------------------------------------------
# Define the organisms and their corresponding strain names and taxonomic IDs.

# A named vector linking organism codes to strain names.
Organisms <- c("73" = "Streptococcus salivariusS73wt",
               "78" = "Streptococcus vestibularisS78wt",
               "93" = "Streptococcus mitisS93wt")

# Extract numeric identifiers from strain names for further use.
Strains <- gsub("[^0-9]", "", Organisms)

# A named vector linking organism codes to NCBI Taxonomy IDs.
taxids <- c("73" = 1304, "78" = 1343, "93" = 28037)

# Create OrgDB (Organism Database) -----------------------------------------
obo_dir <- file.path(here(), "Data", "GO_Check")
go_obo_file <- file.path(obo_dir, "go-basic.obo")
do_not_annotate_file <- file.path(obo_dir, "gocheck_do_not_annotate.obo")

# Create directory if it doesn't exist
if (!dir.exists(obo_dir)) {
  dir.create(obo_dir, recursive = TRUE)
}

# Set up output directories.
output_Quarto <- file.path(here(), "output", "Quarto_Reports")
if (!dir.exists(output_Quarto)) {
  dir.create(output_Quarto, recursive = TRUE)
}

# Path for organism-specific databases.
OrgDBPath <- file.path(here(), "OrgDbs")
if (!dir.exists(OrgDBPath)) {
  dir.create(OrgDBPath, recursive = TRUE)
}



# Obsolete and Non-Annotatable GO Terms -----------------------------------

# Download "gocheck_do_not_annotate.obo" if it's missing
if (!file.exists(do_not_annotate_file)) {
  download.file("https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.obo", 
                destfile = do_not_annotate_file)
}

# Download "go-basic.obo" if it's missing
if (!file.exists(go_obo_file)) {
  download.file("https://purl.obolibrary.org/obo/go/go-basic.obo", 
                destfile = go_obo_file)
}



# Load or Create Organism Databases ----------------------------------------
# Check if organism databases (OrgDB) already exist, and load them. If not, create new ones.

# List of existing organism databases in OrgDBPath.
orgdbpaths <- list.files(OrgDBPath, pattern = "org.*.db", full.names = FALSE)

# If all the required organism databases exist, load them.
if (all(c("org.SmitisS93wt.eg.db", "org.SsalivariusS73wt.eg.db", "org.SvestibularisS78wt.eg.db") %in% orgdbpaths)) {
  for (i in orgdbpaths) {
    library(i, character.only = TRUE)
  }
 } else {
  # If not all databases are present, create them using available data.
  
  # Apply the function to each organism and taxid combination to create the OrgDB.
  mapply(function(taxid, Organism) {
    Strain <- paste(gsub("[^0-9]", "", Organism), "wt", sep = "_")  # Generate strain ID.
    
    # Locate corresponding GO annotation file from Pannzer2 results.
    Pannzer <- list.files(file.path(here(), "Data", "GO_Pannzer2"), pattern = Strain, recursive = TRUE, full.names = TRUE)
    Pannzer <- grep("GO_strain", Pannzer, value = TRUE)

    # If no GO annotation file is found, set Pannzer to NA.
    if (length(Pannzer) == 0) {
      Pannzer <- NA
    }
    
    # Read Bakta annotation file for the strain.
    bakta_table <- read_tsv(list.files(file.path(here(),"Data", "Bakta_annotated_files"),pattern = "tsv",recursive = T,full.names = T), skip = 5)
    colnames(bakta_table) <- gsub("#", "", colnames(bakta_table))  # Clean column names.
    colnames(bakta_table) <- gsub(" ", "_", colnames(bakta_table)) # Replace spaces with underscores.
    
    # Create the OrgDB using Bakta and Pannzer2 data.
    createOrgfromBakta(outdir = OrgDBPath, Strain = Strain, Pannzer = Pannzer, Organism = Organism, 
                       Bakta_tsv = bakta_table, GID = "Locus_Tag", EVIDENCE = "IEA", GOID = "GO", install_DB = FALSE, tax_id = taxid)
    
  }, taxids, Organisms)
  
  # Re-scan the OrgDB directory and install the newly created databases.
  orgdbpaths <- list.files(file.path(here(), "OrgDbs"), pattern = "org.*.db", full.names = T)
  install.packages(orgdbpaths, repos = NULL) 
  lapply(orgdbpaths, require, character.only = TRUE)
}

# Differential Expression Analysis -----------------------------------------
# This section generates differential expression reports using Quarto.

# Define a custom function to render Quarto documents, allowing flexibility in output paths.
process_quarto_render <- function(input, output_dir = NULL, output_file = NULL, ...) {
  x <- quarto::quarto_inspect(input)  # Inspect the Quarto input file.
  output_format <- names(x$formats)  # Extract output formats (e.g., HTML, PDF).
  output_file <- if (is.null(output_file)) x$formats[[output_format]]$pandoc$`output-file` else output_file
  input_dir <- dirname(input)
  output_dir <- if (is.null(output_dir)) input_dir else output_dir
  
  # Render the Quarto document with specified output and parameters.
  quarto::quarto_render(input = input, output_file = output_file, ... = ...)
  
  # Manage file renaming and relocation if necessary.
  output_from <- file.path(input_dir, output_file)
  output_to <- file.path(output_dir, output_file)
  
  if (input_dir != output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (file.copy(output_from, output_to, overwrite = TRUE)) file.remove(output_from)
  } else if (basename(output_from) != basename(output_to)) {
    file.rename(output_from, output_to)
  }
}

# Render differential expression reports for each strain.
for (i in Strains) {
  process_quarto_render("DifferentialExpression.qmd",
                        execute_params = list("strain" = i,"nfcore_output_Path"=nfcore_output_Path),  # Pass the strain as a parameter.
                        output_file = paste0("DifferentialExpression_Strain", i),
                        output_dir = output_Quarto, execute_dir = here())
}

# Cross-Strain Comparison --------------------------------------------------
# Generate a cross-strain comparison report.

process_quarto_render("CrossStrain_Analyis.qmd",
                      output_file = paste0("Cross-strain_Analysis"),
                      output_dir = output_Quarto, execute_dir = here())
