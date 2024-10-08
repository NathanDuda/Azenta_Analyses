

#selected_group_1 <- c('ES.hTSC.CM.')
#selected_group_2 <- c('LS.hTSC.CM.')



DGE_function <- function(selected_group_1, selected_group_2) {
  #aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  aws_prefix <- 'C:/Users/17735/Downloads/Azenta_Analyses/'
  
  
  raw_counts <-  read.csv(paste0(aws_prefix, "Data/Raw_Counts_ProteinCodingGenes.tsv"), sep="") %>%
    select(ID, everything())
  
  
  group_names <- read.table(paste0(aws_prefix, "Data/group_names.txt"), quote="\"", comment.char="")
  replicate_names <- read.table(paste0(aws_prefix, "Data/replicate_names.txt"), quote="\"", comment.char="")
  
  n_replicates <- nrow(replicate_names) / nrow(group_names)
  
  
  g1_group_number <- which(group_names$V1 == selected_group_1)
  g1_col_start <- (1 + (n_replicates * (g1_group_number - 1))) + 1
  g1_col_end <-  (n_replicates * g1_group_number) + 1
  
  g1 <- group_names[g1_group_number, 'V1']
  
  g2_group_number <- which(group_names$V1 == selected_group_2)
  g2_col_start <- (1 + (n_replicates * (g2_group_number - 1))) + 1
  g2_col_end <-  (n_replicates * g2_group_number) + 1
  
  g2 <- group_names[g2_group_number, 'V1']
  
  
  
  g1_colnames <- colnames(raw_counts[g1_col_start:g1_col_end])
  g2_colnames <- colnames(raw_counts[g2_col_start:g2_col_end])
  
  
  
  
  
  
  
  
  # run DESeq2
  rownames(raw_counts) <- NULL
  countData <- raw_counts %>%
    select(all_of(c('ID', g1_colnames, g2_colnames))) %>%
    column_to_rownames(var ='ID')
  
  
  colData <- data.frame(row.names = colnames(countData),
                        type = c(rep(g1, n_replicates), rep(g2, n_replicates)))
  
  
  
  
  library(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ type)
  dds <- DESeq(dds)
  res <- results(dds)
  
  
  
  # make plots
  plot_DGE_function(res, dds, g1, g2)
  
  
  
}




plot_DGE_function <- function(res, dds, g1, g2) {
  #aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  aws_prefix <- 'C:/Users/17735/Downloads/Azenta_Analyses/'
  
  res_df <- as.data.frame(res)
  
  # Add a column to highlight significant genes
  res_df$significance <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, paste0(g2, "_biased"),
                                ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, paste0(g1, "_biased"), "Not Significant"))
  
  res_df <- na.omit(res_df)
  
  
  # MA plot
  #ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  #  geom_point(alpha = 0.5) +
  #  scale_x_log10() +
  #  theme_minimal() +
  #  labs(title = "MA Plot",
  #       x = "Mean of normalized counts",
  #       y = "Log2 Fold Change") +
  #  scale_color_manual(values = c('blue', "red", "gray"))
  
  
  
  
  volcano_plot <- 
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10(p-value)") +
    scale_color_manual(values = c('blue', "red", "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
  
  ggsave(volcano_plot, file = paste0(aws_prefix, 'Data/DGE_volcano_plot.png'), height = 5, width = 7)
  
  
  # Select significantly differentially expressed genes
  sig_genes <- res_df %>%
    filter(significance != 'Not Significant')
  
  # Extract the normalized count data for these genes
  sig_gene_names <- rownames(sig_genes)
  normalized_counts <- assay(rlog(dds))[sig_gene_names, ]
  
  
  row_dist <- dist(normalized_counts)
  row_clust <- hclust(row_dist)
  ordered_rows <- row_clust$order
  
  # Perform hierarchical clustering on the columns (samples)
  col_dist <- dist(t(normalized_counts))
  col_clust <- hclust(col_dist)
  ordered_cols <- col_clust$order
  
  # Reorder the matrix based on the clustering
  clustered_data <- normalized_counts[ordered_rows, ordered_cols]
  
  # Convert the matrix to a long-format data frame
  long_data <- melt(clustered_data)
  
  # Rename columns for clarity
  colnames(long_data) <- c("Gene", "Sample", "Expression")
  
  # Generate the heatmap using ggplot2
  heatmap <- 
    ggplot(long_data, aes(x = Sample, y = Gene, fill = Expression)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("darkblue", "white", 'darkred')) +
    theme_minimal() +
    labs(title = "Heatmap of Significantly DE Genes",
         x = "Sample",
         y = "Gene") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(heatmap, file = paste0(aws_prefix, 'Data/DGE_heatmap.png'), height = 6, width = 5)
  
}


