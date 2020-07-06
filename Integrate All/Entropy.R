library(coop)


#args <-R.utils::commandArgs(asValues=TRUE)
#if (is.null(args[["input"]])) {
#  print("Provide a valid input file name (Batch corrected object) --> RDS file")
#}
#if (is.null(args[["output"]])) {
#  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
#}
#if (is.null(args[["output_matrix"]])) {
#  print("Provide a valid output file name for the matrix. RDS file")
#}



# entropy function
shannon_entropy <- function(x, sample_vector, N_samples) {
  #browser()
  #print(length(sample_vector[x == 1]))
  freq_batch = table(sample_vector[x ==1])/length(sample_vector[x == 1])
  #print(freq_batch)
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_samples))
}

compute_entropy <- function(corrected_space, k_num, dim_num, bool, x, 
                            sample_vector, N_samples,
                            cell_type_vector, N_cell_types, 
                            kit_vector, N_kits, BBKNN=0){
  print(k_num)
  print(dim_num)
  #browser()
  if(BBKNN==1){
    knn_graph <- corrected_space
    knn_graph <- as.matrix(knn_graph)
  } else {
    knn_graph = buildKNNGraph(corrected_space, k=k_num, d=dim_num, transposed = F)[]
    knn_graph <- as.matrix(knn_graph)
  }
  #folder = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Harmony/AllSamples/Batch_Sample_Kit//Cluster/PCA30/res1.4/Entropy/Harmony'
  #write.csv(knn_graph, file = paste0(folder,'_knn_graph.csv'), sep = ",", row.names = T)
  #browser()
  sample_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, sample_vector, N_samples)})
  print('Done sample_entropy')
  celltype_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
  print('Done celltype_entropy')
  kit_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, kit_vector, N_kits)})
  print('Done kit_entropy')
  entropy <- cbind(sample_entropy,celltype_entropy, kit_entropy)
  names(entropy) <- c("sample_entropy", "Cell_type_entropy", "Kit_entropy")
  return(entropy)
}

save_results <- function(x, col_names,folder){
  write.csv(x, file = paste0(folder,'entropy.csv'), sep = ",", row.names = FALSE, col.names = col_names)
}

#2) Seurat objects
compute_entropy_Seurat = function(object, BBKNN = F, BBKNN_KNN = NA,corrected_assay,folder, k_num,dim_num = 10){
  #browser()
  suppressPackageStartupMessages(require(Seurat))
  print("The input object is a Seurat class object")
  
  sample_vector <- object@meta.data[['sample']]
  N_samples <- length(unique(sample_vector))
  
  cell_type_vector <- object$GeneralCellType # Idents(object) #object@meta.data[[celltype_key]]
  N_cell_types <- length(unique(cell_type_vector))
  
  kit_vector <- object@meta.data[['kit']]
  N_kits <- length(unique(kit_vector))
  
  print(unique(sample_vector))
  print(unique(kit_vector))
  print(unique(cell_type_vector))
  
  #write.csv(kit_vector, file = paste0(folder,'_kit_vector.csv'), sep = ",", row.names = T)
  
  #write.csv(cell_type_vector, file = paste0(folder,'_cell_type_vector.csv'), sep = ",", row.names = T)
  
  
  # PCA
  if (!BBKNN){
    tmp = object@reductions[[corrected_assay]]
    tmp = tmp@cell.embeddings
    space <- as.matrix(tmp)
    space = t(space)
  }else{
    space = BBKNN_KNN
    space[space>0] = 1
  }
  folder = paste0(folder,'_k',k_num,'_')
  col_names <- c("sample_entropy", "Cell_type_entropy", "Kit_entropy")
  entropy = compute_entropy(corrected_space = space, 
                            k_num = k_num,dim_num = dim_num, 
                            bool = F, x, 
                            sample_vector = sample_vector, N_samples = N_samples, 
                            cell_type_vector = cell_type_vector,N_cell_types = N_cell_types,
                            kit_vector = kit_vector, N_kits = N_kits , BBKNN = BBKNN)
  save_results(entropy,col_names,folder) 
  print("Entropy calculated over Seurat object!")
}

plotEntropy = function(entropy,folder_output,k_num){
  entropy_melt = melt(entropy)
  pathName <- paste0(folder_output,'_k',k_num,'_Boxplot_entropy','.png')
  png(file=pathName,width=500, height=500, res = 100)
  title = paste0('Median sample: ',specify_decimal(median(entropy$sample_entropy),2), 
                 '\n','Median CellType: ',specify_decimal(median(entropy$celltype_entropy),2), 
                 '\n','Median Kit: ',specify_decimal(median(entropy$kit_entropy),2))
  plot = ggplot(entropy_melt, aes(x=variable, y=value)) + geom_boxplot()
  plot = plot + ggtitle(title) + 
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5))
  print(plot)
  dev.off()
  
  
  pathName <- paste0(folder_output,'_k',k_num,'_Hist_celltype_entropy','.png')
  png(file=pathName,width=500, height=500, res = 100)
  title = paste0('Median: ',specify_decimal(mean(entropy$celltype_entropy), 2) )
  plot = ggplot(entropy, aes(x=celltype_entropy)) + geom_histogram()
  plot = plot + ggtitle(title) + 
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5))
  print(plot)
  dev.off()
  
  pathName <- paste0(folder_output,'_k',k_num,'_Hist_kit_entropy','.png')
  png(file=pathName,width=500, height=500, res = 100)
  title = paste0('Median: ',specify_decimal(mean(entropy$kit_entropy), 2) )
  plot = ggplot(entropy, aes(x=kit_entropy)) + geom_histogram()
  plot = plot + ggtitle(title) + 
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5))
  print(plot)
  dev.off()
  
  pathName <- paste0(folder_output,'_k',k_num,'_Hist_sample_entropy','.png')
  png(file=pathName,width=500, height=500, res = 100)
  title = paste0('Median: ',specify_decimal(mean(entropy$sample_entropy), 2) )
  plot = ggplot(entropy, aes(x=sample_entropy)) + geom_histogram()
  plot = plot + ggtitle(title) + 
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5))
  print(plot)
  dev.off()
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))