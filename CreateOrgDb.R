# Load Required Libraries --------------------------------------------------
library(here)
source(file.path(here(), "Libraries.R"))


# Create DB from Bakta ----------------------------------------------------
exchangeCOGValues <- function(df) {
  # Ensure the dataframe has the columns we need
  if (!("COG" %in% names(df)) || !("COG_category" %in% names(df))) {
    stop("The dataframe does not have the required columns (COG, COG_category).")
  }

  # Find rows where COG_category contains "COG"
  rowsToExchange <- grepl("COG", df$COG_category)

  # Exchange values for those rows
  temp <- df$COG[rowsToExchange]
  df$COG[rowsToExchange] <- df$COG_category[rowsToExchange]
  df$COG_category[rowsToExchange] <- temp

  # Return the modified dataframe
  return(df)
}


createOrgfromBakta <- function(Bakta_tsv, outdir, Organism, Strain, Pannzer = NA, GID = "query", EVIDENCE = "IEA", GOID = "GOs", install_DB = T, tax_id) {
  # Convert comma-separated DbXrefs to long format
  data_long <- Bakta_tsv %>%
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
    mutate(Gene = ifelse(!is.na(Gene), paste0(Gene), paste0(Product))) %>%
    dplyr::select(dplyr::all_of(GID), everything())

  # Clean up the data frame
  data_wide$"NA" <- NULL
  data_wide <- data_wide[complete.cases(data_wide[GID]), ]

  # Remove columns that contain only NA values
  allna <- apply(data_wide, 2, function(x) all(is.na(x)))
  data_wide <- data_wide[, !allna]

  # Remove specific database columns
  data_wide <- data_wide[, !colnames(data_wide) %in% c("BlastRules", "IS", "VFDB", "RFAM")]

  # Format GO terms with GO: prefix
  data_wide$GO <- ifelse(is.na(data_wide$GO), NA, paste0("GO:", gsub("; ", ", GO:", data_wide$GO)))
  Bakta_tsv <- data_wide

  # Prepare column names
  columnstoinclude <- colnames(Bakta_tsv)
  columnstoinclude <- columnstoinclude[!columnstoinclude == GID]
  colnames(Bakta_tsv)[colnames(Bakta_tsv) == GID] <- "GID"

  # Create gene information data frame
  gene_info <- data.frame(
    GID = Bakta_tsv$GID,
    SYMBOL = Bakta_tsv$Gene,
    DESCRIPTION = Bakta_tsv$Product
  )

  # Process GO terms
  goterms <- gsub(" ", "", as.character(Bakta_tsv$GO))
  go_terms <- strsplit(goterms, ",")
  gene2go <- data.frame(
    GID = rep(Bakta_tsv$GID, sapply(go_terms, length)),
    GO = unlist(go_terms),
    EVIDENCE = EVIDENCE
  )

  # Process Pannzer annotations if provided
  if (!is.na(Pannzer)) {
    # Helper function to replace last digits in GO terms
    replace_last_digits <- function(go, num) {
      n <- nchar(num)
      substr(go, 1, nchar(go) - n) <- substr(go, 1, nchar(go) - n)
      paste0(substr(go, 1, nchar(go) - n), num)
    }

    # Read and process Pannzer annotations
    GOannotations <- read.csv(Pannzer, sep = "\t")
    GOannotations$GO <- paste("GO:0000000")
    GOannotations$GO <- replace_last_digits(GOannotations$GO, GOannotations$goid)
    GOannotations$EVIDENCE <- EVIDENCE

    # Combine Pannzer annotations with existing GO terms
    Pannzer_GO <- GOannotations[c("qpid", "GO", "EVIDENCE")]
    colnames(Pannzer_GO) <- c("GID", "GO", "EVIDENCE")
    gene2go <- rbind(gene2go, Pannzer_GO)
    gene2go <- gene2go %>% distinct()
  }

  # Load and check GO ontology
  obo_dir <- file.path(here(), "Data", "GO_Check")
  go_obo_file <- file.path(obo_dir, "go-basic.obo")
  do_not_annotate_file <- file.path(obo_dir, "gocheck_do_not_annotate.obo")

  # Filter out obsolete and do-not-annotate terms
  do_not_annotate <- get_ontology(do_not_annotate_file)
  goBasic <- get_ontology(go_obo_file)
  obsolete <- goBasic$id[goBasic$obsolete]
  gene2go <- gene2go[!gene2go$GO %in% c(do_not_annotate$id, obsolete), ]

  # Create KEGG annotation data frames
  gene2ko <- data.frame(
    GID = Bakta_tsv$GID,
    KEGG = Bakta_tsv$KEGG
  )
  
  # Create COG annotation data frames
  gene2cog <- data.frame(
    GID = Bakta_tsv$GID,
    COGID = Bakta_tsv$COG,
    COGCAT = Bakta_tsv$COG_category
  )

  # Remove empty entries from annotation data frames
  gene2go <- gene2go[!is.na(gene2go$GO) & gene2go$GO != "", ]
  gene2cog <- gene2cog[!is.na(gene2cog$COGID) & gene2cog$COGID != "", ]
  gene2ko <- gene2ko[!is.na(gene2ko$KEGG) & gene2ko$KEGG != "", ]

  # Create the organism package using AnnotationForge
  My_OrgDB <- AnnotationForge::makeOrgPackage(
    gene_info = gene_info,
    goTable = "go",
    go = gene2go,
    ko = gene2ko,
    cog = gene2cog,
    maintainer = "Bernd Daller <Bernd.Daller@ukr.de>",
    author = "Bernd Daller",
    outputDir = outdir,
    tax_id = tax_id,
    genus = str_split(Organism, pattern = " ")[[1]][1],
    species = str_split(Organism, pattern = " ")[[1]][2],
    verbose = T,
    version = "0.01"
  )

  # Install the created database if requested
  if (install_DB) {
    orgdbpath <- list.files(outdir, pattern = "org.*.db", full.names = T)
    install.packages(orgdbpath, repos = NULL)
  }
}
