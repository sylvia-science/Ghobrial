PipelineIntegrateAll = function(data_integrate,sample_name,folder_output,sampleParam,integrate_merge){
  file_str = ''
  #browser()
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)
  
  
  
  ########################
  ########################
  
  if (integrate_merge == 'Merge'){
    #browser()
    #####################
    ## QC
    #####################
    #nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
    #                          ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
    #percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
    #data = quality_control(data,filter,nFeature_RNA_list,percent_mt,sample_name)
    data_integrate = NormalizeData(data_integrate, normalization.method = "LogNormalize", scale.factor = 10000)
    
    ########################
    # Get Variable Genes
    ########################
    
    nfeatures_val = 2000 #sampleParam$nfeatures_val[sampleParam['Sample'] == sample_name]
    data_integrate = FindVariableFeatures(data_integrate, selection.method = "vst", nfeatures = nfeatures_val)
    
    dir.create( paste0(folder_output,'QC Metrics/'), recursive = TRUE)
    
    pathName = paste0(folder_output,'QC Metrics/FindVariableFeatures','.png')
    print(pathName)
    png(file=pathName,width=600, height=350, res = 100)
    print(VariableFeaturePlot(data_integrate) + ylim(0,10))
    dev.off()
    

  }
  
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  
  cell_features = getCellMarkers()
  
  #Score for cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data_integrate = CellCycleScoring(data_integrate, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  data_integrate = ScaleData(data_integrate, features = rownames(data_integrate))
  
  print(paste0('PCA: ', PCA_dim))
  data_integrate <- RunPCA(data_integrate, features = VariableFeatures(object =data_integrate), npcs = PCA_dim)
  
  ######################
  ## Start Plotting
  ######################
  
  #Visualize PCA results
  visualize_PCA(data_integrate,folder_output,PCA_dim)
  

  pathName = paste0(folder_output,'PCA/elbow',file_str,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data_integrate,ndims = PCA_dim))
  dev.off()
  
  # Jackstraw
  #data_integrate <- JackStraw(data_integrate, num.replicate = 100)
  #data_integrate <- ScoreJackStraw(data_integrate, dims = 1:PCA_dim)
  
  
  # Find Neighbors and Clusters
  data_integrate <- getCluster (data_integrate,resolution_val, PCA_dim)
  
  # Find Cluster Biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data_integrate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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
  print(DimPlot(data_integrate, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data_integrate, features = top10$gene))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  print(FeaturePlot(data_integrate, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_split', file_str,'.png'))  
  png(file=pathName,width=600, height=350)
  print(DimPlot(data_integrate, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
  dev.off()
  
  return(data_integrate)
  
}