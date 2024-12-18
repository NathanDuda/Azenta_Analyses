

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
#gene_lengths <- read.delim("C:/Users/17735/Downloads/Azenta_Analyses/Azenta_Projects/30-1005983745/hit-counts/ES-CM1.counts.txt")


#project_path <- list(path = './Azenta_Projects/30-1041694139/', project_type = 'azenta')
#project_path <- list(path = './Azenta_Projects/30-1005983745/', project_type = 'azenta')
#project_path <- list(path = './Genomics_Projects/HS24900/', project_type = 'genomics')

normalization_function <- function(project_path, normalization_type, exp_cutoff, final_id_format = 'NA') {
  #aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  aws_prefix <- 'C:/Users/17735/Downloads/Azenta_Analyses/'
  

  if (project_path$project_type == 'azenta') {
    raw_counts <- read.csv(file.path(project_path$path, 'hit-counts', 'raw_counts.csv'))
    #############raw_counts <- read.csv('./Azenta_Projects/30-1041694139/hit-counts/raw_counts.csv')
    # select gene lengths file
    file_list <- list.files(paste0(project_path$path, '/hit-counts/'))
    file_list <- file_list[!file_list %in% c('raw_counts.csv', 'TPM_values.csv')]
    gene_length_path <- paste0(project_path$path, '/hit-counts/', file_list[1])
    
    # import the gene lengths file
    gene_lengths <- read.csv(gene_length_path, sep = "")
    
    
    ENSG_symbol <- raw_counts %>%
      select(ID, Symbol = Gene.name)
    
    raw_counts <- raw_counts %>% select(-Gene.name)
    
  }
  
  if (project_path$project_type == 'genomics') {
    # import raw_counts file
    raw_counts <- read.delim(file.path(project_path$path, 'Fulltable.txt'))
    ############raw_counts <- read.delim('./Genomics_Projects/HS24900/Fulltable.txt')
    
    # format raw_counts file to match azenta
    raw_counts <- raw_counts %>% select(-starts_with("TPM_"))
    colnames(raw_counts)[1] <- 'ID'
    
    # select gene lengths file
    file_list <- list.files(paste0(project_path$path, '/'))
    file_list <- file_list[!file_list %in% c('Fulltable.txt')]
    gene_length_path <- paste0(project_path$path, '/', file_list[1])
    gene_lengths <- read.delim(gene_length_path, comment.char="#")
    
    colnames(gene_lengths)[1] <- 'Geneid'
  }
  
  
  
  # import pre-existing data 
  protein_coding_ENSGs <- read.csv(paste0(aws_prefix, "PreExisting_Data/protein_coding_ENSG_ids.tsv"), sep = "")
  
  # keep only protein coding genes
  raw_counts <- raw_counts %>%
    filter(ID %in% protein_coding_ENSGs$ENSG) %>%
    select(ID, everything())
  
  # write protein coding genes to file
  write.table(raw_counts, file = paste0(aws_prefix, 'Data/Raw_Counts_ProteinCodingGenes.tsv'))
  
  
  # format the input files 
  gene_lengths <- gene_lengths %>%
    select(ID = Geneid, Length)
  
  counts <- raw_counts %>%
    merge(., gene_lengths, by = 'ID') %>%
    select(ID, Length, everything())
  
  # get the names of each group of replicates
  colnames(counts) <- gsub("\\.[^.]*\\.$", "", colnames(counts)) # () turned into .. in the column names, if a colname ends in . (the genomics project does) then remove the characters in between so that the colname will end in the replicate 
  group_names <- colnames(counts)
  replicate_names <- group_names[!group_names %in% c('ID','Length')]
  
  
  # make it work for the triple replicates labeled 1,2,nothing instead of 1,2,3
  if (all(grepl("(1|2|3)$", replicate_names))) {
    # all values end in 1, 2, or 3  -> remove last character 
    group_names <- lapply(replicate_names, function(x) substr(x, 1, nchar(x) - 1))
  } else {
    # Use sub to remove only the last dot and numbers that follow (only if only numbers follow)
    group_names <- sub("\\.[0-9]+$", "", replicate_names)
  }  
  
  
  # get list of group names 
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
  counts <- counts %>% select(ID, Length, all_of(ordered_replicate_names$V1))
  
  norm_counts <- counts 
  for (group in colnames(counts)[3:ncol(counts)]) {
    if (normalization_type == 'RPKM') {
      total_reads <- sum(counts[[group]])
      total_reads <- total_reads / 1e6 
      
      norm_counts <- norm_counts %>%
        mutate(length_kb = Length / 1000,
               !!paste0('norm_', group) := .data[[group]] / (length_kb * total_reads))
    }
    if (normalization_type == 'TPM') {
      norm_counts <- norm_counts %>% mutate(n_over_l = .data[[group]] / Length)
      denominator <- sum(norm_counts$n_over_l) / 1e6
      
      norm_counts <- norm_counts %>%
        mutate(!!paste0('norm_', group) := n_over_l / denominator)
    }
  }
  
  # keep only normalized columns
  norm_counts <- norm_counts %>%
    select(ID, starts_with(paste('norm_')))
  
  # write normalized counts to file
  write.table(norm_counts, file = paste0(aws_prefix, 'Data/', normalization_type, '_normalized_exp.tsv'))
  
  
  
  # get average expression of replicates in each group 
  for (i in 1:length(group_names$V1)) {
    
    # get group name - becomes name of averaged column 
    name <- group_names$V1[i]
    
    # determine the starting and ending columns to average 
    start_col <- (i - 1) * n_replicates + 2
    end_col <- start_col + n_replicates - 1
    
    # calculate the average 
    norm_counts <- norm_counts %>%
      mutate(!!sym(name) := rowMeans(across(all_of(start_col:end_col))))
    
  }
  
  # select only the averaged columns
  averaged_exp <- norm_counts %>%
    select(all_of(c('ID', group_names$V1)))
  
  # after averaging, set expression values lower than exp_cutoff to 0 
  averaged_exp <- averaged_exp %>%
    mutate(across(2:ncol(averaged_exp), ~ if_else(. <= exp_cutoff, 0, .)))
  
  
  # write averaged expression to file
  write.table(averaged_exp, file = paste0(aws_prefix, 'Data/Averaged_', normalization_type, '_exp.tsv'))
  
  
  if (final_id_format == 'symbol') {
    symbol_averaged_exp <- averaged_exp %>%
      merge(., ENSG_symbol, by = 'ID') %>%
      select(Symbol, ID, everything())
    
    return(symbol_averaged_exp)
  }
  
  if (final_id_format == 'ENSG') {
    return(averaged_exp)
  }
  
}

  