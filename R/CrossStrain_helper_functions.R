#' Helper functions for cross-strain analysis in RNAseq data
#' 
#' This file contains reusable functions for differential gene expression analysis

#' Process COG categories based on regulation direction
#'
#' This function extracts and processes COG categories from differential expression data
#' based on regulation direction (upregulated or downregulated) and significance threshold.
#'
#' @param DE_table Data frame containing differential expression results
#' @param log2FC_threshold Threshold for log2 fold change filtering
#' @param direction Direction of regulation, either "up" or "down"
#' @return Data frame with COG categories and their counts
#' @examples
#' # upregulated <- process_cog_categories(de_table, 0.5, "up")
process_cog_categories <- function(DE_table, log2FC_threshold, direction = "up") {
  # Set the comparison operator based on direction
  comp_operator <- if(direction == "up") `>` else `<`
  threshold <- if(direction == "up") log2FC_threshold else -log2FC_threshold
  
  # Extract and process COG categories
  DE_table<-DE_table$COG_category[DE_table$padj < 0.05 & 
                       comp_operator(DE_table$log2FoldChange, threshold) & 
                       complete.cases(DE_table$COG_category)] %>%
    trimws() %>%
    str_split("") %>%
    unlist() %>%
    table() %>%
    as.data.frame()
  colnames(DE_table) <- c("COG_Category", "Count")
  return(DE_table)
}

#' Generate radar plots for COG data
#'
#' Creates a radar plot for visualizing COG category distributions by strain/condition group
#'
#' @param data Data frame containing COG categories and their values
#' @param group_name String indicating which group to subset from the data
#' @param object_colors Named vector of colors for each group
#' @return A plot_grid object containing the radar plot and its legend
#' @examples
#' # radar_plots <- generate_radar_plot(cog_data, "S. salivarius", strain_colors)
generate_radar_plot <- function(data, group_name, object_colors) {
  group_colours <- object_colors[grepl(group_name, names(object_colors))]

  # Subset the data for the specific group
  ggradardata <- subset(data, grepl(group_name, data$group))
  
  # Generate the radar plot for upregulated genes
  radar_regulated <- ggradardata[ggradardata$Regulation == "Upregulated", -1] %>%
    ggradar(
      grid.label.size = 2,
      axis.label.size = 2,
      group.point.size = 2,
      legend.text.size = 5,
      grid.min = 0,
      grid.mid = max(ggradardata[, -c(1:2)]) / 6,
      grid.max = max(ggradardata[, -c(1:2)]) / 3,
      label.gridline.max = "25%",
      gridline.mid.colour = "black",
      label.gridline.min = "0%",
      label.gridline.mid = "10%",
      background.circle.colour = "#1A120B",
      group.line.width = 1,
      group.colours = group_colours,
      fill = TRUE, fill.alpha = 0.3, legend.position = "right"
    )
    
  # Extract and arrange legend with plot
  legend <- get_legend(radar_regulated +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 2, lineheight = 0.9, face = "bold"),
      legend.box = "horizontal",
      legend.background = element_rect(fill = "#FAF9F6", linewidth = 0.5, linetype = "blank")
    ))
    
  plot_grid(radar_regulated + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")), 
            legend, 
            ncol = 2)
}

#' Extract differentially expressed genes based on threshold
#'
#' This function extracts gene names from a gene expression dataset that meet
#' specific fold change thresholds.
#'
#' @param gene_data Named vector of log fold changes
#' @param threshold Log2 fold change threshold
#' @param direction Direction of regulation ("up" or "down")
#' @return Vector of gene names that meet the threshold criteria
#' @examples
#' # upregulated_genes <- extract_de_genes(gene_set, 0.5, "up")
extract_de_genes <- function(gene_data, threshold, direction = "up") {
  if(direction == "up") {
    filtered <- gene_data[gene_data > threshold]
  } else {
    filtered <- gene_data[gene_data < -threshold]
  }
  
  out <- unique(names(filtered))
  out <- out[!is.na(out)]
  return(out)
}

#' Generate UpSet plot for gene sets
#'
#' Creates a UpSet plot 
#'
#' @param gene_matrix Matrix from make_comb_mat for the gene sets
#' @param sample_colors Named vector of colors for each sample
#' @param pt_size Size of the points in the plot (in mm)
#' @param font_size Font size for labels
#' @return An UpSet plot object
#' @examples
#' # upset_plot <- create_upset_plot(upregulated_matrix, sample_colors, 3, 5)
create_upset_plot <- function(gene_matrix, sample_colors, pt_size = 3, font_size = 5) {
  ComplexHeatmap::UpSet(gene_matrix,
    set_order = names(sample_colors),
    row_names_gp = gpar(fontsize = font_size),
    column_names_gp = gpar(fontsize = font_size),
    lwd = 0.5,
    comb_order = order(comb_size(gene_matrix), decreasing = TRUE),
    pt_size = unit(pt_size, "mm"),
    bg_pt_col = sample_colors,
    bg_col = rev(sample_colors),
    comb_col = c("#b69978", "#927552", "#72583d", "#5a3f31", "#2d2722", "black")[comb_degree(gene_matrix)],
    top_annotation = upset_top_annotation(gene_matrix, 
                                         height = unit(3, "cm"), 
                                         axis_param = list(gp = gpar(fontsize = font_size)), 
                                         annotation_name_gp = gpar(fontsize = font_size + 1)),
    right_annotation = upset_right_annotation(gene_matrix,
      width = unit(1.5, "cm"),
      gp = gpar(fill = sample_colors),
      annotation_name_gp = gpar(fontsize = font_size + 1),
      annotation_name_side = "top",
      axis_param = list(side = "top", gp = gpar(fontsize = font_size))
    )
  )
}

#' Process set data from Venn diagrams
#'
#' Processes data from Venn diagram regions into a tidy format
#'
#' @param venn_data Data from a Venn diagram object
#' @return A data frame with region details and items
#' @examples
#' # overlap_data <- process_region_data(Venn(upregulated_genes))
process_set_data <- function(venn_data) {
  result <- data.frame(
    region = names(venn_data@sets),
    count = sapply(venn_data@sets, length),
    item = sapply(venn_data@sets, function(x) paste(x, collapse = ", "))
  )
  return(result)
}
