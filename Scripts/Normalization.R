
# assumes every group has the same number of replicates 
# make sure tpm and rpkm values are actually accurate 

library(tidyverse)

# indicate the type of normalization
normalization_type <- 'RPKM'

# provide the raw_counts file found in the hit-counts folder
raw_counts <- read.csv("C:/Users/17735/Downloads/raw_counts.csv")

# provide the expression file for any sample found in the hit-counts folder
# this is only for the gene lengths, so it doesn't matter which one
gene_lengths <- read.delim("C:/Users/17735/Downloads/ES-CM1.counts.txt")


# import pre-existing data 
protein_coding_ENSGs <- read.csv("./PreExisting_Data/protein_coding_ENSG_ids.tsv", sep="")



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
replicate_factor <- factor(replicate_groups, levels = group_names)
ordered_replicate_names <- replicate_names[order(replicate_factor, replicate_names)]

# check that each group has the same number of 
t <- table(replicate_groups)
if(min(t) != max(t)) {paste0('ERROR: The number of replicates varies accross groups.')}

# get the number of replicates
n_replicates <- length(replicate_names) / length(group_names)

# reorder the counts dataframe so that every n_replicates columns is a new group
counts <- counts %>% select(ID, Length, Gene.name, all_of(ordered_replicate_names))

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


# keep only protein coding genes
norm_counts <- norm_counts %>%
  filter(ID %in% protein_coding_ENSGs$ENSG)

# write normalized counts to file
write.table(norm_counts, file = paste0('./Data/', normalization_type, '_normalized_counts.tsv'))













