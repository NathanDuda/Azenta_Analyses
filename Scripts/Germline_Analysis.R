

germline_function <- function(normalization_type){
  aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  
  germline_data <- read.csv(paste0(aws_prefix, "PreExisting_Data/germline_data.tsv"), sep="") %>%
    select(ID = ENSG, everything())
  
  exp <- read.csv(paste0(paste0(aws_prefix, 'Data/Averaged_'), normalization_type, '_exp.tsv'), sep="")
  
  group_names <- read.table(paste0(aws_prefix, "Data/group_names.txt"), quote="\"", comment.char="")
  
  
  # get results 
  germline_results <- exp %>%
    merge(., germline_data, by = 'ID') %>% 
    filter(rowSums(select(., 3:(nrow(group_names) + 2)) == 0) != 4)
  
  # if output exists, write to file 
  if (nrow(germline_results) > 0) {
    write.table(germline_results, file = paste0(aws_prefix, 'Data/Germline_Results.tsv'))
    
    # plot results
    plot_upset_germline_function(germline_results = germline_results, group_names = group_names)
  }
  if (nrow(germline_results) == 0) {
    write.table('No expressed genes were found in any germline', file = paste0(aws_prefix, 'Data/Germline_Results.tsv'))
  }
  
  germline_table <- germline_results %>%
    pivot_longer(cols = (ncol(germline_results) -2):ncol(germline_results), names_to = "column", values_to = "value") %>%
    group_by(ID) %>%
    filter(value == 1) %>%
    summarize(combined = paste(column, collapse = ", ")) %>%
    merge(germline_results, by = "ID") %>%
    select(ID, Gene.name, germline_layer = combined)
  
  return(germline_table)
}


plot_upset_germline_function <- function(germline_results, group_names){
  aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
  
  categories <- c('mesoderm', 'endoderm', 'neuroectoderm')
  
  
  binary_matrix <- germline_results %>%
    pivot_longer(cols = 3:(nrow(group_names) + 2)) %>%
    filter(value != 0) %>%
    select(neuroectoderm, endoderm, mesoderm) %>%
    as.data.frame()
  
  
  
  ##binary_matrix <- sapply(categories, function(x) grepl(x, germline_results$germline))
  #binary_matrix <- as.data.frame(t(as.data.frame(binary_matrix)))
  #colnames(binary_matrix) <- categories
  
  #binary_matrix <- t(binary_matrix)
  #rownames(binary_matrix) <- germline_results$ID
  #binary_matrix[binary_matrix == T] <- 1
  

  upset_matrix <- make_comb_mat(binary_matrix)
  

  # save plot to image
  #png("./Data/Germline_upset_plot.png", width = 4.5, height = 3, units = 'in', res = 1200)
  #ComplexHeatmap::UpSet(upset_matrix)
  #dev.off()
  
  
  
  library(ggplotify)
  library(fs)
  
  ggplot2::ggsave(ggplotify::as.ggplot(UpSet(upset_matrix, comb_order = rev(order(comb_size(upset_matrix))))),
                  filename = fs::path(paste0(aws_prefix, "Data/Germline_upset_plot.png")),
                  device = "png",
                  units = "in",
                  height = 2.5, width = 4.5)
  
}

















