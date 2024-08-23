

pathway_function <- function(normalization_type){
  
  germline_data <- read.csv("./PreExisting_Data/germline_data.tsv", sep="") %>%
    select(ID = ENSG, Gene.name = gn_symbol, germline)
  
  exp <- read.csv(paste0('./Data/Averaged_', normalization_type, '_exp.tsv'), sep="")
  
  group_names <- read.table("./Data/group_names.txt", quote="\"", comment.char="")
  
  
  # get results 
  germline_results <- exp %>%
    merge(., germline_data, by = c('ID', 'Gene.name')) %>% 
    filter(rowSums(select(., 2:(nrow(group_names) + 2)) == 0) != 4)
  
  # if output exists, write to file 
  if (nrow(pathway_results > 0)) {
    write.table(germline_results, file = './Data/Germline_Results.tsv')
  }
  if (nrow(pathway_results == 0)) {
    print('No expressed genes were found in any pathway')
  }
}


