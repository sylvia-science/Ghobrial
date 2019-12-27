RunPipelineIntegrateAll = function(data_integrated,sample_name,folder_output,sampleParam){
  file_str = ''
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)
  
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  
  ########################
  ########################
  cell_features = getCellMarkers()
  
  
  #Score for cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data_integrated = CellCycleScoring(data_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  data_integrated = ScaleData(data_integrated, features = rownames(data_integrated))
  
  print(paste0('PCA: ', PCA_dim))
  data_integrated <- RunPCA(data_integrated, features = VariableFeatures(object =data_integrated), npcs = PCA_dim)
  
  ######################
  ## Start Plotting
  ######################
  
  #Visualize PCA results
  visualize_PCA(data_integrated,folder_output,PCA_dim)
  

  pathName = paste0(folder_output,'PCA/elbow',file_str,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data_integrated,ndims = PCA_dim))
  dev.off()
  
  # Jackstraw
  #data_integrated <- JackStraw(data_integrated, num.replicate = 100)
  #data_integrated <- ScoreJackStraw(data_integrated, dims = 1:PCA_dim)
  
  
  # Find Neighbors and Clusters
  data_integrated <- getCluster (data_integrated,resolution_val, PCA_dim)
  
  # Find Cluster Biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  all_markers =  markers %>% group_by(cluster)
  
  
  
  all_markers$Cell = NA
  # Add known markers to all
  #browser()
  
  for (i in 1:nrow(all_markers)){
    marker_rows = grep(all_markers$gene[i], cell_features$Markers, value=TRUE)
    marker_idx = which( cell_features$Markers %in% marker_rows)
    #browser()
    if (length(marker_idx) > 0){
      #browser()
      cell_list = cell_features$Cell[marker_idx]
      cell_list = paste(cell_list, sep="", collapse=", ") 
      all_markers$Cell[i] = cell_list
    }
  }

  write.csv(all_markers, file = paste0(folder_output,'AllFeatures',file_str,'.csv'),row.names=FALSE)
  
  
  # Visualize clustering
  filepath_cluster = paste0( folder_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data_integrated, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data_integrated, features = top10$gene))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  print(FeaturePlot(data_integrated, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_split', file_str,'.png'))  
  png(file=pathName,width=600, height=350)
  print(DimPlot(data_integrated, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
  dev.off()
  
  return(data_integrated)
  
}