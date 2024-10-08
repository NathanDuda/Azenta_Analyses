

pathway_function <- function(normalization_type){
  #aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  aws_prefix <- 'C:/Users/17735/Downloads/Azenta_Analyses/'
  
  pathway_data <- read.csv(paste0(aws_prefix, "PreExisting_Data/pathway_data.tsv"), sep="") %>%
    select(ID = ENSG, pathway)
  
  exp <- read.csv(paste0(aws_prefix, 'Data/Averaged_', normalization_type, '_exp.tsv'), sep="")
  
  group_names <- read.table(paste0(aws_prefix, "Data/group_names.txt"), quote="\"", comment.char="")
  
  
  # get results 
  pathway_results_original <- exp %>%
    merge(., pathway_data, by = c('ID')) %>% 
    filter(rowSums(select(., 2:(nrow(group_names) + 2)) == 0) != 4)
  
  # if output exists, write to file 
  if (nrow(pathway_results_original) > 0) {
    
    pathway_results <- pathway_results_original %>%
      select(-ID) %>%
      pivot_longer(cols = 1:(ncol(.) - 1)) %>%
      filter(value != 0) %>%
      select(pathway, name)
    
    pathway_results <- as.data.frame(table(pathway_results$name, pathway_results$pathway)) %>%
      filter(Freq != 0) %>%
      select(line = Var1, pathway = Var2, n_genes = Freq)
    
    write.table(pathway_results, file = paste0(aws_prefix, 'Data/Pathway_Results.tsv'))
  }
  if (nrow(pathway_results) == 0) {
    write.table('No expressed genes were found in any pathway', file = paste0(aws_prefix, 'Data/Pathway_Results.tsv'))
  }
  
  pathway_table <- pathway_results_original %>%
    pivot_longer(cols = 3:(ncol(.) - 1)) %>%
    filter(value != 0) %>%
    select(Line = name, ID, pathway) %>%
    distinct()
  
  return(pathway_table)
}

