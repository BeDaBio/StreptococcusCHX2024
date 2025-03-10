#' Exchange COG values and categories
#'
#' This function exchanges the COG values with COG categories when the COG_category
#' column contains "COG" entries, correcting a common issue in COG annotations.
#'
#' @param df A data frame containing COG and COG_category columns
#'
#' @return A data frame with corrected COG and COG_category values
#' @export
#'
#' @examples
#' # df <- exchange_cog_values(gene_infos)
exchange_cog_values <- function(df) {
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
  
  return(df)
}

#' Process DESeq2 Results
#'
#' This function processes the results from a DESeq2 analysis, creates output directories,
#' formats the results, and saves them as Excel files.
#'
#' @param res A DESeqResults object from DESeq2
#' @param comparison Character string naming the comparison being analyzed
#' @param results_path Path where results will be saved
#' @param postfix Optional string to append to output filenames
#' @param gene_infos Optional data frame with gene annotation information
#' @param by Vector of length 2 specifying which columns to use for merging gene_infos
#'
#' @return A data frame containing significant results (padj < 0.05)
#' @export
#'
#' @examples
#' # process_deseq(res_CHXresponse_WT, results_path = results_path,
#' #               comparison = "CHX125.vs.Untreated_inWTcells", 
#' #               gene_infos = gene_infos, by = c("gene", "Locus Tag"))
process_deseq <- function(res, 
                          comparison, 
                          results_path, 
                          postfix = "", 
                          gene_infos = NULL, 
                          by = c("gene", "reference")) {
  
  # Load required packages if not already loaded
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required but not installed")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed")
  }
  
  # Create directory for comparison if it doesn't exist
  if (!dir.exists(file.path(results_path, comparison))) {
    dir.create(file.path(results_path, comparison), recursive = TRUE)
  }

  # Convert DESeq results to a data frame and arrange by padj
  all <- res %>%
    as.data.frame() %>%
    arrange(padj)
  all$gene <- gsub("\\..*", "", rownames(all))
  all$comparison <- comparison

  # Optionally merge with gene information
  if (!length(gene_infos) == 0) {
    all <- merge(
      as.data.frame(all), 
      gene_infos, 
      by.x = by[1], 
      by.y = by[2], 
      sort = FALSE, 
      all.x = TRUE, 
      all.y = FALSE
    )
  }

  # Save full DESeq results
  openxlsx::write.xlsx(
    all, 
    file.path(results_path, comparison, paste0(comparison, "_Differential_Expression_all", postfix, ".xlsx")), 
    rowNames = FALSE
  )

  # Filter significant results (padj < 0.05)
  fdr <- res %>%
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(desc(log2FoldChange))

  if (nrow(fdr) == 0) {
    openxlsx::write.xlsx(
      fdr, 
      file.path(results_path, comparison, paste0(comparison, "_Differential_Expression_fdr0.05", postfix, ".xlsx")), 
      rowNames = FALSE
    )
    return(fdr)
  } else {
    fdr$gene <- gsub("\\..*", "", rownames(fdr))
    fdr$comparison <- comparison
    if (!length(gene_infos) == 0) {
      fdr <- merge(
        as.data.frame(fdr), 
        gene_infos, 
        by.x = by[1], 
        by.y = by[2], 
        sort = FALSE, 
        all.x = TRUE, 
        all.y = FALSE
      )
    }
    openxlsx::write.xlsx(
      fdr, 
      file.path(results_path, comparison, paste0(comparison, "_Differential_Expression_fdr0.05", postfix, ".xlsx")), 
      rowNames = FALSE
    )
    return(fdr)
  }
}

#' Create a Volcano Plot from Differential Expression Results
#'
#' This function reads differential expression results from an Excel file
#' and creates an enhanced volcano plot highlighting the top differentially
#' expressed genes.
#'
#' @param x Path to the Excel file containing differential expression results
#' @param remove_pattern Regex pattern to remove from the comparison name (default: "_.*")
#'
#' @return An EnhancedVolcano plot object
#' @export
#'
#' @examples
#' # volcano_plot <- plot_volcano("path/to/differential_expression_results.xlsx")
plot_volcano <- function(x, remove_pattern = "_.*") {
  # Load required packages if not already loaded
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required but not installed")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed")
  }
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    stop("Package 'EnhancedVolcano' is required but not installed")
  }
  
  y <- read_xlsx(x)
  
  if (length(colnames(y)) < 2) {
    print("nothing happened here")
    return(NA)
  }
  
  # Extract name, title, and comparison from filename
  name <- gsub("\\.xlsx", "", basename(x))
  title <- gsub("_Differential_Expression.*", "", name)
  comparison <- gsub(remove_pattern, "", name)
  comparisons <- strsplit(comparison, "\\.vs\\.")[[1]]
  
  y <- y[, c("gene_id", "pvalue", "log2FoldChange", "pvalue")]
  y <- y[complete.cases(y), ]
  
  # Select top and bottom 20 genes by log2FoldChange
  selected_genes <- y[y$log2FoldChange < (-1), ]
  bottom_genes <- selected_genes[order(selected_genes$pvalue), ] %>% head(10)
  selected_genes <- y[y$log2FoldChange > 1, ]
  top_genes <- selected_genes[order(selected_genes$pvalue), ] %>% head(10)
  
  #avoid duplicated gene_ids
  y[["gene_id"]][duplicated(y[["gene_id"]])] <- paste(
    y[["gene_id"]][duplicated(y[["gene_id"]])],
    seq_along(y[["gene_id"]][duplicated(y[["gene_id"]])])
  )
  
  # Generate Volcano Plot
  plot <- EnhancedVolcano(
    y,
    lab = y[["gene_id"]],
    x = "log2FoldChange",
    y = "pvalue",
    title = NULL,
    axisLabSize = 10,
    labSize = 2,
    selectLab = c(top_genes$gene_id, bottom_genes$gene_id),
    max.overlaps = Inf,
    legendLabels = c(
      "NS",
      expression(Log[2] ~ FC ~ " > 2"),
      "FDR < 0.05",
      expression("FDR < 0.05" ~ and ~ log[2] ~ FC)
    ),
    caption = NULL,
    subtitle = NULL,
    drawConnectors = TRUE
  )
  
  
  
  return(plot)

}

#' Plot normalized gene expression values
#'
#' This function creates plots of normalized gene expression values from differential 
#' expression analysis results, allowing for visualization of specific genes across 
#' different conditions.
#'
#' @param de_table Path to Excel file containing differential expression results
#' @param res Results object from DESeq2
#' @param comparison Character string naming the comparison being analyzed
#' @param int_group Vector of factors to use for grouping in plotCounts
#' @param x_value Column to use on the x-axis
#' @param facet Column to use for faceting
#' @param selected_genes Vector of gene IDs to plot
#' @param dds DESeq2 dataset object
#' @param metadata Data frame containing sample metadata
#' @param genotype_chx_colors Named vector of colors for different genotype/CHX combinations
#' @param theme_custom Custom theme function for plot styling (optional)
#'
#' @return A nested list of ggplot2 objects
#' @export
#'
#' @examples
#' # plot_list <- plot_normalized_values(de_table = "path/to/results.xlsx", 
#' #                                    res = res_object,
#' #                                    comparison = "comparison_name",
#' #                                    int_group = c("Genotype", "CHX"),
#' #                                    x_value = "CHX",
#' #                                    facet = "Genotype",
#' #                                    selected_genes = c("gene1", "gene2"),
#' #                                    dds = dds_object,
#' #                                    metadata = metadata,
#' #                                    genotype_chx_colors = color_vector)
plot_normalized_values <- function(de_table, 
                                  res, 
                                  comparison, 
                                  int_group = c("Treatment", "Celltype"), 
                                  x_value, 
                                  facet, 
                                  selected_genes,
                                  dds,
                                  metadata,
                                  genotype_chx_colors,
                                  theme_custom = NULL) {
  # Load required packages if not already loaded
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required but not installed")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed")
  }
  
  
  # Read data from Excel file
  de_table <- readxl::read_xlsx(de_table)
  if (length(colnames(de_table)) < 2) {
    print("nothing happened here")
    return(NA)
  }
  facet <- rlang::ensym(facet)

  # Process each selected gene
  sapply(selected_genes, function(genetype) {
    title <- genetype 
    de_table <- de_table %>%
      dplyr::arrange(desc(log2FoldChange)) %>%
      dplyr::filter(`#Sequence Id` == "contig_1")
    
    locus <- de_table[(grepl(genetype, de_table$"gene", fixed = TRUE) | 
                     grepl(genetype, de_table[, "gene_id"], fixed = TRUE)) & 
                     de_table$log2FoldChange > 0.05, "gene"]
    
    title <- de_table[(grepl(genetype, de_table$"gene", fixed = TRUE) | 
                     grepl(genetype, de_table[, "gene_id"], fixed = TRUE)) & 
                     de_table$log2FoldChange > 0.05, "gene_id"]
                     
    wrapped_title <- sapply(title, function(x) paste(strwrap(x, width = 50), collapse = "\n"))

    # Process each locus for the current gene
    sapply(locus, function(x) {
      FoldChange <- format(as.numeric(de_table[de_table$"gene" == x, "log2FoldChange"]), digits = 2)
      Padj <- format(as.numeric(de_table[de_table$"gene" == x, "padj"]), digits = 2)
      
      # Get normalized counts for the current gene
      Counts <- DESeq2::plotCounts(dds, x, intgroup = int_group, returnData = TRUE, transform = TRUE)
      
      # Merge with metadata
      Counts <- merge(
        Counts, 
        metadata[, !colnames(metadata) %in% colnames(Counts)], 
        by.x = "row.names", 
        by.y = "SRR", 
        all.x = TRUE, 
        all.y = FALSE
      )

      # Generate plot
      p <- ggplot2::ggplot(
        Counts, 
        aes(
          x = CHX, 
          y = count, 
          shape = Genotype, 
          color = Genotype_CHX, 
          fill = Genotype_CHX
        )
      ) +
        ggplot2::stat_summary(
          fun.y = mean, 
          geom = "crossbar", 
          size = 0.5, 
          aes(group = Genotype_CHX)
        ) + 
        ggplot2::geom_jitter(
          color = "black",
          size = 3,
          position = ggplot2::position_dodge2(width = 0.3)
        ) +
        ggplot2::facet_wrap(~Genotype) +
        ggplot2::theme_bw() + 
        theme_custom() +
        ggplot2::scale_shape_manual(values = c(21, 24)) +
        ggplot2::scale_fill_manual(values = genotype_chx_colors) +
        ggplot2::scale_color_manual(values = genotype_chx_colors) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            order = 1,
            title.position = "left", 
            title.hjust = 0,
            legend.title.align = 0,
            override.aes = list(
              size = 3, 
              linetype = c(0, 0, 0, 0), 
              shape = c(17, 17, 16, 16), 
              color = c(genotype_chx_colors[order(names(genotype_chx_colors))]), 
              fill = c(NA, NA, NA, NA)
            )
          ),
          shape = "none", 
          color = "none"
        ) +
        ggplot2::labs(
          title = wrapped_title,
          x = NULL,
          y = "counts"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 7)) + 
        ggplot2::expand_limits(y = c(0, ceiling(max(Counts$count) + 0.2 * max(Counts$count)))) 
        
      return(p)
    }, simplify = FALSE)
  }, simplify = FALSE)
}

#' Perform Gene Set Enrichment Analysis (GSEA) with GO terms
#'
#' This function performs Gene Set Enrichment Analysis using Gene Ontology terms
#' from a differential expression results file and handles filtering of obsolete terms.
#'
#' @param input_path Path to Excel file containing differential expression results
#' @param output_path Path where results will be saved (default: same directory as input)
#' @param ont Ontology to use ("BP", "MF", or "CC")
#' @param org_db OrgDb object for the organism
#'
#' @return A simplified GSEA results object
#' @export
#'
#' @examples
#' # gse_results <- go_enriched_clusterprofiler(input_path = "path/to/results.xlsx",
#' #                                           ont = "BP",
#' #                                           org_db = org.myOrganism.db)
go_enriched_clusterprofiler <- function(input_path, 
                                       output_path = NULL, 
                                       ont, 
                                       org_db) {
  # Load required packages if not already loaded
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required but not installed")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required but not installed")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required but not installed")
  }
  if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
    stop("Package 'ontologyIndex' is required but not installed")
  }

  name <- gsub("\\.xlsx", "", basename(input_path))
  de_table <- readxl::read_xlsx(input_path)
  de_table <- de_table[!duplicated(de_table["gene_id"]), ]

  if (length(colnames(de_table)) < 1) {
    print("nothing happened here")
    return(NA)
  }

  de_table <- de_table %>% dplyr::arrange(desc(.data[["log2FoldChange"]]))
  genes <- de_table %>% dplyr::pull(.data[["log2FoldChange"]], "gene_id")

  genes <- genes[!is.na(names(genes)) & !is.na(genes)]
  if ("comparison" %in% colnames(de_table)) {
    comparison <- unique(de_table$comparison)
  } else {
    comparison <- ""
  }

  sample_name <- gsub("\\.xlsx", "", basename(input_path))

  # Check input and set up output directories
  if (is.null(output_path)) {
    output_path <- dirname(input_path)
  }

  if (!dir.exists(file.path(output_path, "GO_GSEA", ont))) {
    dir.create(
      file.path(output_path, paste0("GO_GSEA_", ont)), 
      recursive = TRUE, 
      showWarnings = FALSE
    )
  }

  # Perform GSEA
  gse <- tryCatch(
    {
      gse <- clusterProfiler::gseGO(
        geneList = genes,
        ont = ont,
        keyType = "SYMBOL",
        minGSSize = 2,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE,
        OrgDb = org_db,
        pAdjustMethod = "BH",
        seed = TRUE
      )
      saveRDS(gse, file = file.path(output_path, paste0(sample_name, "GO_GSEA.rds")))
      gse
    },
    error = function(cond) {
      warning(paste0("Error occurred during gseGO: ", cond$message))
      return(NA)
    }
  )
  
  # Process results if successful
  if(inherits(gse, "gseaResult")) {
    # Get paths for GO ontology files
    # Allow for customization via environment variables
    go_check_dir <- file.path(here(), "Data", "GO_Check")
    print(go_check_dir)
    # Load the ontology files
    tryCatch({
      do_not_annotate <- ontologyIndex::get_ontology(file.path(go_check_dir, "gocheck_do_not_annotate.obo"))
      goBasic <- ontologyIndex::get_ontology(file.path(go_check_dir, "go-basic.obo"))
      obsolete <- goBasic$id[goBasic$obsolete]
      
      # Filter out obsolete terms
      fine <- gse@result[!gse@result$ID %in% c(obsolete, do_not_annotate$id), ]
      gse@result <- fine
      
      # Simplify and save results
      gse2 <- clusterProfiler::simplify(gse)
      openxlsx::write.xlsx(gse, file = file.path(output_path, paste0(sample_name, "_GSEA.xlsx")))
      openxlsx::write.xlsx(gse2, file = file.path(output_path, paste0(sample_name, "_GSEA_nonRedundant.xlsx")))
      return(gse2)
    }, error = function(e) {
      warning("Could not filter obsolete GO terms: ", e$message)
      # Save unsimplified results if ontology files couldn't be loaded
      openxlsx::write.xlsx(gse, file = file.path(output_path, paste0(sample_name, "_GSEA.xlsx")))
      return(gse)
    })
  } else {
    return(gse)
  }
}
