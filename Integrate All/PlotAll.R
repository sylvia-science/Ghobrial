plotAll = function(data,folder_base, folder,sample_name,sampleParam,label_TF,
                   integrate_TF = FALSE, # This is an old variable, keep it can stay false
                   DE_perm_TF = FALSE, # If you want to do DE between cluster and each individual cluster
                   clusterTF = F, # If you want to recluster (Good for if you are changing the resolution val)
                   markersTF = T, # If you want to do DE between each cluster and all other clusters combined.
                   keepOldLabels = F, # If you want to force the program to keep old labels
                   groupBy = NA, # Varibles that will be used to make umap plots where the variables are split in the plot
                   PCA_dim = NA, # If NA, will gram PCA_dim from sampleParam
                   resolution_val = NA # If NA, will gram resolution_val from sampleParam
                   ){
  #source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
  require(gtools)
  browser()
  print(folder)
  print(paste0('label_TF: ', label_TF))
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
  if (is.na(PCA_dim)){
    PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  }
  if (is.na(resolution_val)){
    resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  }
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
  # Score for cell cycle
  s.genes = cc.genes$s.genes
  g2m.genes = cc.genes$g2m.genes
  data = CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  # Known Cell Markers
  cell_features = getCellMarkers(folder_base)
  
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
  #visualize_PCA(data,folder,PCA_dim)

  dir.create( paste0(folder,'PCA'), recursive = TRUE)
  pathName <- paste0(folder,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()

  #############
  ## TEST: REFIND VAR FEATURE
  #############

  #browser()
  nfeatures_val = 2000 
  #data_RNA = data
  #DefaultAssay(object = data_RNA) <- "RNA"
  #data_RNA@active.assay = 'RNA'
  #data_post_varFeatures = FindVariableFeatures(data_RNA, selection.method = "vst", nfeatures = nfeatures_val)
  
  #browser()
  
  # Cluster with Umap
  if (clusterTF == TRUE){
    data = getCluster (data,resolution_val, PCA_dim)
  }else if (keepOldLabels == F){
    #data = FindNeighbors(data, dims = 1:PCA_dim)
    data = FindClusters(data, resolution = resolution_val)
  }
  #browser()
  
  # Name cells
  if (label_TF & keepOldLabels == F){
    #browser()
    cluster_IDs <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
    data = label_cells(data,cluster_IDs)
    # Percentage of cell type
  }
  
  ##############################
  ## Get Cluster Stats
  ##############################
  
  cluster_num = clusterStats(data)
  print(cluster_num)
  
  filepath_stats = paste0(filepath_cluster, 'Stats/')
  dir.create( filepath_stats, recursive = TRUE)
  write.csv(cluster_num, file = paste0(filepath_stats,'clusterStats',file_str,'.csv'),row.names = FALSE)
  
  
  if (any(which(data$orig.ident == "data_pre"))){
    data_pre = data[which(data$orig.ident == "data_pre")]
    cluster_num_pre = clusterStats(data_pre)
    write.csv(cluster_num_pre, file = paste0(filepath_stats,'clusterStats_pre',file_str,'.csv'),row.names = FALSE)
    
  }
  
  if (any(which(data$orig.ident == "data_post"))){
    data_post = data[which(data$orig.ident == "data_post")]
    cluster_num_post = clusterStats(data_post)
    write.csv(cluster_num_post, file = paste0(filepath_stats,'clusterStats_post',file_str,'.csv'),row.names = FALSE)
    
  }
  
  ########################
  ## Visualize clustering
  ########################
  
  
  pathName <- paste0(filepath_cluster,
                     paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1000, res = 100)
  print(DimPlot(data,pt.size = 0.5, reduction = "umap",label = TRUE))
  dev.off()
  

  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data,pt.size = 0.5, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'percent.mt','.png')
  png(file=pathName,width=600, height=600)
  print(FeaturePlot(data,pt.size = 0.5, features = c("percent.mt"), split.by = "split_var"))
  dev.off()
  #browser()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap','_PCA',PCA_dim,resolution_val,'_splitAll', '','.png'))  
  png(file=pathName,width=2600, height=500,res = 100)
  print(DimPlot(data, label=T, repel=F, reduction = "umap", split.by = "split_var"))
  dev.off()
  
  if (integrate_TF){
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_split',file_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
    dev.off()
  }
  
  if (!is.na(groupBy)){
    for (group in groupBy){
      pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy',group,'.png'))
      png(file=pathName,width=1000, height=1000)
      plot = DimPlot(data,pt.size = 0.5, reduction = "umap",label = FALSE,group.by  = group)
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20)
        
      )
      print(plot)
      
      dev.off()
      
      pathName <- paste0(filepath_cluster,paste0('ClusterUmap','_PCA',PCA_dim,'_res',resolution_val,'_split', group,'.png'))  
      png(file=pathName,width=3000, height=1000)
      plot = DimPlot(data, label=T, repel=F, reduction = "umap", split.by = group)
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20)
        
      )
      print(plot)
      dev.off()
      
    }
  }
  browser()
  
  ##################################
  
  ##############################
  ## Find Cluster Biomarkers
  ##############################
  if (markersTF == TRUE){
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    num_markers = 10
    markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers  %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
    # Plotting the top 10 markers for each cluster.
    top10 = markers %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
    all_markers =  markers %>% group_by(cluster)
    
    pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
    png(file=pathName,width=1000, height=1200)
    print(DoHeatmap(data, features = top10$gene))
    dev.off()
    
    
    # Add known markers
    all_markers = cellMarkers(all_markers,cell_features)
    #write.csv(all_markers, file = paste0(filepath_cluster,'AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
    write.csv(all_markers, file = paste0(filepath_cluster,'Features',file_str,'.csv'),row.names=FALSE)
  } 
  
  # Get DE between cluster permutations
  if (DE_perm_TF){
    filepath_mvn = paste0( filepath_cluster, 'DE/nVsm/')
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
                        ,'.csv')
      print(filepath)
      #browser()
      write.csv(df_output, file = filepath,row.names = FALSE)
    }
  }

}


