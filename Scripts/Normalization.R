
# assumes every group has the same number of replicates 
# make sure tpm and rpkm values are actually accurate 

#library(tidyverse)

# indicate the type of normalization
#normalization_type <- 'RPKM'

# indicate the cutoff expression value
#exp_cutoff <- 1

# provide the raw_counts file found in the hit-counts folder
#raw_counts <- read.csv("C:/Users/17735/Downloads/raw_counts.csv")

# provide the expression file for any sample found in the hit-counts folder
# this is only for the gene lengths, so it doesn't matter which one
#gene_lengths <- read.delim("C:/Users/17735/Downloads/ES-CM1.counts.txt")



normalization_function <- function(raw_counts, gene_lengths, normalization_type, exp_cutoff) {
  aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  
  # import pre-existing data 
  protein_coding_ENSGs <- read.csv(paste0(aws_prefix, "PreExisting_Data/protein_coding_ENSG_ids.tsv"), sep="")
  
  
  # keep only protein coding genes
  raw_counts <- raw_counts %>%
    filter(ID %in% protein_coding_ENSGs$ENSG) %>%
    select(ID, Gene.name, everything())
  
  # write protein coding genes to file
  write.table(raw_counts, file = paste0(aws_prefix, 'Data/Raw_Counts_ProteinCodingGenes.tsv'))
  
  
  # format the input files 
  gene_lengths <- gene_lengths %>%
    select(ID = Geneid, Length)
  
  counts <- raw_counts %>%
    merge(., gene_lengths, by = 'ID') %>%
    select(ID, Length, Gene.name, everything())
  
  # get the names of each group of replicates
  group_names <- colnames(counts)
  replicate_names <- group_names[!group_names %in% c('ID','Length','Gene.name')]
  group_names <- lapply(replicate_names, function(x) substr(x, 1, nchar(x) - 1))
  group_names <- as.character(unique(group_names))
  
  # order the replicate names based on their groups 
  replicate_groups <- sapply(strsplit(replicate_names, "\\d+"), `[`, 1)
  write.table(as.data.frame(replicate_groups), file = paste0(aws_prefix, 'Data/replicate_groups.txt')) # only for unit test
  replicate_factor <- factor(replicate_groups, levels = group_names)
  ordered_replicate_names <- replicate_names[order(replicate_factor, replicate_names)]
  
  # write group list and replicate group list to files
  ordered_replicate_names <- as.data.frame(ordered_replicate_names)
  colnames(ordered_replicate_names) <- 'V1'
  write.table(ordered_replicate_names, file = paste0(aws_prefix, "Data/replicate_names.txt"))

  group_names <- as.data.frame(group_names)
  colnames(group_names) <- 'V1'
  write.table(group_names, file = paste0(aws_prefix, "Data/group_names.txt"))
  
  # check that each group has the same number of replicates
  t <- table(replicate_groups)
  if(min(t) != max(t)) {#stop("The number of replicates varies accross groups.")
    return()}
  
  # get the number of replicates
  n_replicates <- nrow(ordered_replicate_names) / nrow(group_names)
  
  # reorder the counts dataframe so that every n_replicates columns is a new group
  counts <- counts %>% select(ID, Length, Gene.name, all_of(ordered_replicate_names$V1))
  
  # TPM normalize if chosen 
  if (normalization_type == 'TPM') {
    norm_counts <- counts %>%
      mutate(across(4:(ncol(counts)), ~ .x / (Length / 1000), .names = "_{col}")) %>%
      mutate(across(starts_with("_"), ~ .x / sum(.x) * 1e6, .names = "norm{col}")) %>%
      select(ID, Gene.name, starts_with("norm_"))
  }
  
  # RPKM normalize if chosen 
  if (normalization_type == 'RPKM') {
    norm_counts <- counts %>%
      mutate(across(4:(ncol(counts)), ~ .x / (Length / 1000), .names = "_{col}")) %>%
      mutate(across(starts_with("_"), ~ .x / (sum(get(gsub("_", "", cur_column()))) / 1e6), .names = "norm{col}")) %>%
      select(ID, Gene.name, starts_with("norm_"))
  }
  
  # set expression values lower than exp_cutoff to 0 
  exp <- norm_counts %>%
    mutate(across(3:ncol(norm_counts), ~ if_else(. <= exp_cutoff, 0, .)))
  
  # write normalized counts to file
  write.table(exp, file = paste0(aws_prefix, 'Data/', normalization_type, '_normalized_exp.tsv'))
  
  
  
  # get average expression of replicates in each group 
  for (i in 1:length(group_names$V1)) {
    
    # get group name - becomes name of averaged column 
    name <- group_names$V1[i]
    
    # determine the starting and ending columns to average 
    start_col <- (i - 1) * n_replicates + 3
    end_col <- start_col + n_replicates - 1
    
    # calculate the average 
    exp <- exp %>%
      mutate(!!sym(name) := rowMeans(across(all_of(start_col:end_col))))
    
  }
  
  # select only the averaged columns
  averaged_exp <- exp %>%
    select(all_of(c('ID', 'Gene.name', group_names$V1)))
  
  # after averaging, set expression values lower than exp_cutoff to 0 
  averaged_exp <- averaged_exp %>%
    mutate(across(3:ncol(averaged_exp), ~ if_else(. <= exp_cutoff, 0, .)))
  
  
  # write averaged expression to file
  write.table(averaged_exp, file = paste0(aws_prefix, 'Data/Averaged_', normalization_type, '_exp.tsv'))
  
  
}
