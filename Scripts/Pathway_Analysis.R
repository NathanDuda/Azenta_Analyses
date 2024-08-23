

pathway_function <- function(normalization_type){
  
  pathway_data <- read.csv("./PreExisting_Data/pathway_data.tsv", sep="") %>%
    select(ID = ENSG, Gene.name = gn_symbol, pathway)
  
  exp <- read.csv(paste0('./Data/Averaged_', normalization_type, '_exp.tsv'), sep="")
  
  group_names <- read.table("./Data/group_names.txt", quote="\"", comment.char="")
  
  
  # get results 
  pathway_results <- exp %>%
    merge(., pathway_data, by = c('ID', 'Gene.name')) %>% 
    filter(rowSums(select(., 2:(nrow(group_names) + 2)) == 0) != 4)
  
  # if output exists, write to file 
  if (nrow(pathway_results) > 0) {
    
    pathway_results <- pathway_results %>%
      select(-ID, -Gene.name) %>%
      pivot_longer(cols = 1:(ncol(.) - 1)) %>%
      filter(value != 0) %>%
      select(pathway, name)
    
    pathway_results <- as.data.frame(table(pathway_results$name, pathway_results$pathway)) %>%
      filter(Freq != 0) %>%
      select(line = Var1, pathway = Var2, n_genes = Freq)
    
    write.table(pathway_results, file = './Data/Pathway_Results.tsv')
  }
  if (nrow(pathway_results) == 0) {
    write.table('No expressed genes were found in any pathway', file = './Data/Pathway_Results.tsv')
  }
  return(pathway_results)
}









