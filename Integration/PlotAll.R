plotAll = function(data,folder,sample_name,sampleParam,label_TF,integrate_TF = FALSE,DE_perm_TF = FALSE,clusterTF = TRUE){
  source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
  require(gtools)
  #browser()
  print(folder)
  print(paste0('label_TF: ', label_TF))
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
  
  # Known Cell Markers
  cell_features = getCellMarkers()
  
  #cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  if (label_TF){
    file_str = '_label'
    
  }else
  {
    file_str = ''
  }
  print(paste0('file_str: ', file_str))
  
  ##################
  ## PCA
  ##################
  #browser()
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  visualize_PCA(data,folder,PCA_dim)

  pathName <- paste0(folder,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  ########################
  # Plot Variable Genes
  ########################
  # pathName = paste0(folder,'QC Metrics/FindVariableFeatures','.png')
  # print(pathName)
  # png(file=pathName,width=600, height=350, res = 100)
  # print(VariableFeaturePlot(data) + ylim(0,10))
  # dev.off()
  
  # Cluster with Umap
  ## ADD DOTPLOT
  resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  #clusterTF = TRUE
  if (clusterTF == TRUE){
    data = getCluster (data,resolution_val, PCA_dim)
  }
  
  
  # Name cells
  if (label_TF){
    #browser()
    cluster_IDs <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
    #browser()
    data = label_cells(data,cluster_IDs)
    # Percentage of cell type
  }
  
  ##############################
  ## Get Cluster Stats
  ##############################
  
  cluster_num = clusterStats(data@active.ident)
  print(cluster_num)
  write.csv(cluster_num, file = paste0(folder,'Stats/clusterStats',file_str,'.csv'),row.names = FALSE)
  
  
  if (any(which(data$orig.ident == "data_pre"))){
    data_pre = data@active.ident[which(data$orig.ident == "data_pre")]
    cluster_num_pre = clusterStats(data_pre)
    write.csv(cluster_num_pre, file = paste0(folder,'Stats/clusterStats_pre',file_str,'.csv'),row.names = FALSE)
    
  }
  
  if (any(which(data$orig.ident == "data_post"))){
    data_post = data@active.ident[which(data$orig.ident == "data_post")]
    cluster_num_post = clusterStats(data_post)
    write.csv(cluster_num_post, file = paste0(folder,'Stats/clusterStats_post',file_str,'.csv'),row.names = FALSE)
    
  }
  
  ##############################
  ## Find Cluster Biomarkers
  ##############################
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  all_markers =  markers %>% group_by(cluster)
  
  # Add known markers
  all_markers = cellMarkers(all_markers,cell_features)
  
  write.csv(all_markers, file = paste0(folder,'DE/AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  
  # Get DE between cluster permutations
  if (DE_perm_TF){
    filepath_mvn = paste0( folder, 'DE/nVsm/')
    print(filepath_mvn)
    dir.create( filepath_mvn, recursive = TRUE)
    Features_mvn_df = diffExpPermute(data,cell_features)
    
    level_list = unique(levels(data))
    for (level in 1:length(level_list) ){
      ident1 = level_list[level]
      
      if (ident1 == '?'){
        ident1 = 'NA'
      }
      
      df_output = Features_mvn_df[[level]]
      df_output = df_output[c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","ident_1","ident_2","gene", "Cell")]
      filepath = paste0(filepath_mvn
                        ,'Features_',ident1,'Vsn',file_str
                        ,'_Patient',Patient_num,'.csv')
      print(filepath)
      write.csv(df_output, file = filepath,row.names = FALSE)
    }
  }
  ########################
  ## Visualize clustering
  ########################
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
  write.csv(all_markers, file = paste0(filepath_cluster,'AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data, features = top10$gene))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  
  if (integrate_TF){
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_split',file_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
    dev.off()
  }
}


