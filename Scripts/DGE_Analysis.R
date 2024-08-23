


#selected_group_1 <- c('ES.hTSC.CM.')
#selected_group_2 <- c('LS.CM')



DGE_function <- function(selected_group_1, selected_group_2) {
  
  
  raw_counts <-  read.csv("./Data/Raw_Counts_ProteinCodingGenes.tsv", sep="") %>%
    select(ID, Gene.name, everything())
  
  
  group_names <- read.table("./Data/group_names.txt", quote="\"", comment.char="")
  replicate_names <- read.table("./Data/replicate_names.txt", quote="\"", comment.char="")
  
  n_replicates <- nrow(replicate_names) / nrow(group_names)
  
  
  g1_group_number <- which(group_names$V1 == selected_group_1)
  g1_col_start <- (g1_group_number * (n_replicates - (g1_group_number-1))) + 2
  g1_col_end <- (g1_group_number * n_replicates) + 2
  
  g1 <- group_names[g1_group_number, 'V1']
  
  g2_group_number <- which(group_names$V1 == selected_group_2)
  g2_col_start <- (g2_group_number * n_replicates - (g2_group_number-1)) + 2
  g2_col_end <- (g2_group_number * n_replicates) + 2
  
  g2 <- group_names[g2_group_number, 'V1']
  
  
  
  g1_colnames <- colnames(raw_counts[g1_col_start:g1_col_end])
  g2_colnames <- colnames(raw_counts[g2_col_start:g2_col_end])
  
  
  
  
  
  
  
  
  # run DESeq2
  rownames(raw_counts) <- NULL
  countData <- raw_counts %>%
    select(all_of(c('ID', g1_colnames, g2_colnames))) %>%
    column_to_rownames(var ='ID')
  
  
  colData <- data.frame(row.names = colnames(t),
                        type = c(rep(g1, n_replicates), rep(g2, n_replicates)))
  
  
  
  
  library(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ type)
  dds <- DESeq(dds)
  res <- results(dds)
  
  
  
  
  
  
  
}
















