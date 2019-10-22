plotAll = function(data,folder,sample_name){
  print(sample_name)
  print(folder)
  
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  metaData <- read_excel(filename_metaData)
  sampleParam <- read_excel(filename_sampleParam)
  

  #Get Variable Genes
  norm_val <- sampleParam$norm_val[sampleParam['Sample'] == sample_name]
  nfeatures_val <- sampleParam$nfeatures_val[sampleParam['Sample'] == sample_name]
  return_list= gene_var(data,norm_val,nfeatures_val)
  top10 = return_list$top10
  

  # PCA
  visualize_PCA(data,folder,sample_name)
  visualize_dim(data)
  #JackStrawPlot(data, dims = 1:15)
  
  plot = ElbowPlot(data)
  pathName <- paste0(folder,'Cluster/elbow.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  # Cluster with Umap
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  plot = DimPlot(data, reduction = "umap",label = TRUE)
  pathName <- paste0(folder,paste0('Cluster/ClusterUmap',resolution_val,'.png'))
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  # Find Cluster Biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  #Visualize clustering
  plot = DoHeatmap(data, features = top10$gene)
  pathName <- paste0(folder,'Cluster/HeatMap.png')
  png(file=pathName,width=1000, height=1200)
  print(plot)
  dev.off()
  
  #Visualize clustering
  
  pathName <- paste0(folder,'Cluster/ClusterMetrics.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  # Get gene Descriptions
  #gene_desc_top10 =  get_gene_desc(top10)
  #gene_desc_top10
  return(data)
  
}
