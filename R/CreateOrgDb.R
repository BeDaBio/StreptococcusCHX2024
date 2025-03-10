# Load Required Libraries --------------------------------------------------
library(here)
source(file.path(here(), "R/Libraries.R"))


# Create DB from Bakta ----------------------------------------------------

#' Exchange COG and COG_category values
#'
#' @param df A data frame containing COG and COG_category columns
#' @return A data frame with exchanged values for rows where COG_category contains "COG"
#' @examples
#' df <- data.frame(COG = c("A", "B"), COG_category = c("COG123", "X"))
#' exchangeCOGValues(df)
exchangeCOGValues <- function(df) {
  # Check for required columns
  required_cols <- c("COG", "COG_category")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Find rows where COG_category contains "COG"
  rows_to_exchange <- grepl("COG", df$COG_category)

  # Exchange values for those rows
  temp <- df$COG[rows_to_exchange]
  df$COG[rows_to_exchange] <- df$COG_category[rows_to_exchange]
  df$COG_category[rows_to_exchange] <- temp

  return(df)
}


#' Create organism database from Bakta TSV output
#'
#' @param bakta_tsv Data frame from Bakta TSV file
#' @param outdir Directory to save the organism database
#' @param organism Organism name (e.g., "Streptococcus mutans")
#' @param strain Strain name
#' @param pannzer Path to Pannzer annotation file (optional)
#' @param gid_col Column name for gene ID (default: "query")
#' @param evidence_code Evidence code for GO annotations (default: "IEA")
#' @param go_col Column name for GO terms (default: "GOs")
#' @param install_db Whether to install the database (default: TRUE)
#' @param tax_id Taxonomy ID of the organism
#' @return Path to the created organism database
createOrgFromBakta <- function(bakta_tsv, outdir, organism, strain, 
                              pannzer = NA, gid_col = "query", 
                              evidence_code = "IEA", go_col = "GOs", 
                              install_db = TRUE, tax_id) {
  
  # Convert comma-separated DbXrefs to long format
  data_long <- bakta_tsv %>%
    mutate(DbXrefs = strsplit(as.character(DbXrefs), ",\\s*")) %>%
    unnest(DbXrefs)

  # Split DbXrefs into Database and Reference columns
  data_long <- data_long %>%
    mutate(
      Database = sub(":.*", "", DbXrefs),
      Reference = sub(".*:", "", DbXrefs)
    ) %>%
    dplyr::select(-DbXrefs)

  # Group references by other columns and concatenate multiple references
  data_grouped <- data_long %>%
    group_by(across(-Reference), Database) %>%
    summarise(Reference = paste(unique(Reference), collapse = "; "), .groups = "drop")

  # Convert to wide format and process COG annotations
  data_wide <- data_grouped %>%
    pivot_wider(names_from = Database, values_from = Reference) %>%
    separate_wider_delim(COG, delim = ";", names = c("COG", "COG_category")) %>%
    exchangeCOGValues() %>%
    mutate(Gene = ifelse(!is.na(Gene), Gene, Product)) %>%
    dplyr::select(dplyr::all_of(gid_col), everything())

  # Clean up the data frame
  data_wide$"NA" <- NULL
  data_wide <- data_wide[complete.cases(data_wide[gid_col]), ]

  # Remove columns that contain only NA values
  na_cols <- apply(data_wide, 2, function(x) all(is.na(x)))
  data_wide <- data_wide[, !na_cols]

  # Remove specific database columns
  exclude_cols <- c("BlastRules", "IS", "VFDB", "RFAM")
  data_wide <- data_wide[, !colnames(data_wide) %in% exclude_cols]

  # Format GO terms with GO: prefix
  data_wide$GO <- ifelse(is.na(data_wide$GO), NA, paste0("GO:", gsub("; ", ", GO:", data_wide$GO)))
  bakta_processed <- data_wide

  # Prepare column names
  columns_to_include <- colnames(bakta_processed)
  columns_to_include <- columns_to_include[columns_to_include != gid_col]
  colnames(bakta_processed)[colnames(bakta_processed) == gid_col] <- "GID"

  # Create gene information data frame
  gene_info <- data.frame(
    GID = bakta_processed$GID,
    SYMBOL = bakta_processed$Gene,
    DESCRIPTION = bakta_processed$Product
  )

  # Process GO terms
  go_terms_raw <- gsub(" ", "", as.character(bakta_processed$GO))
  go_terms_list <- strsplit(go_terms_raw, ",")
  gene2go <- data.frame(
    GID = rep(bakta_processed$GID, sapply(go_terms_list, length)),
    GO = unlist(go_terms_list),
    EVIDENCE = evidence_code
  )

  # Process Pannzer annotations if provided
  if (!is.na(pannzer)) {
    # Helper function to replace last digits in GO terms
    replace_last_digits <- function(go, num) {
      prefix_length <- nchar(go) - nchar(num)
      paste0(substr(go, 1, prefix_length), num)
    }

    # Read and process Pannzer annotations
    go_annotations <- read.csv(pannzer, sep = "\t")
    go_annotations$GO <- "GO:0000000"
    go_annotations$GO <- mapply(replace_last_digits, 
                               go_annotations$GO, 
                               go_annotations$goid)
    go_annotations$EVIDENCE <- evidence_code

    # Combine Pannzer annotations with existing GO terms
    pannzer_go <- go_annotations[c("qpid", "GO", "EVIDENCE")]
    colnames(pannzer_go) <- c("GID", "GO", "EVIDENCE")
    gene2go <- rbind(gene2go, pannzer_go)
    gene2go <- distinct(gene2go)
  }

  # Load and check GO ontology
  obo_dir <- file.path(here(), "Data", "GO_Check")
  go_obo_file <- file.path(obo_dir, "go-basic.obo")
  do_not_annotate_file <- file.path(obo_dir, "gocheck_do_not_annotate.obo")

  # Filter out obsolete and do-not-annotate terms
  do_not_annotate <- get_ontology(do_not_annotate_file)
  go_basic <- get_ontology(go_obo_file)
  obsolete <- go_basic$id[go_basic$obsolete]
  gene2go <- gene2go[!gene2go$GO %in% c(do_not_annotate$id, obsolete), ]

  # Create KEGG annotation data frames
  gene2ko <- data.frame(
    GID = bakta_processed$GID,
    KEGG = bakta_processed$KEGG
  )
  
  # Create COG annotation data frames
  gene2cog <- data.frame(
    GID = bakta_processed$GID,
    COGID = bakta_processed$COG,
    COGCAT = bakta_processed$COG_category
  )

  # Remove empty entries from annotation data frames
  gene2go <- gene2go[!is.na(gene2go$GO) & gene2go$GO != "", ]
  gene2cog <- gene2cog[!is.na(gene2cog$COGID) & gene2cog$COGID != "", ]
  gene2ko <- gene2ko[!is.na(gene2ko$KEGG) & gene2ko$KEGG != "", ]

  # Create the organism package using AnnotationForge
  genus_species <- str_split(organism, pattern = " ")[[1]][1:2]
  
  my_org_db <- AnnotationForge::makeOrgPackage(
    gene_info = gene_info,
    goTable = "go",
    go = gene2go,
    ko = gene2ko,
    cog = gene2cog,
    maintainer = "Bernd Daller <Bernd.Daller@ukr.de>",
    author = "Bernd Daller",
    outputDir = outdir,
    tax_id = tax_id,
    genus = genus_species[1],
    species = genus_species[2],
    verbose = TRUE,
    version = "0.01"
  )

  # Install the created database if requested
  if (install_db) {
    org_db_path <- list.files(outdir, pattern = "org.*.db", full.names = TRUE)
    if (length(org_db_path) > 0) {
      install.packages(org_db_path, repos = NULL)
      return(org_db_path)
    } else {
      warning("No database package found to install")
    }
  }
  
  return(my_org_db)
}
