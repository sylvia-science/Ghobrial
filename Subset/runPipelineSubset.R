runPipelineSubset = function(data,folder_input,sample_name,sampleParam,filter,scale_TF){
  require(gtools)
  print(sample_name)
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
  
  cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
  #browser()
  file_str = ''
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)

  ########################
  ## Known Cell Markers
  ########################
  cell_features = getCellMarkers()
  
  ############
  ## QC
  # Do we need to normalize for subclustering?
  ############
  #data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  ########################
  # Get Variable Genes
  ########################
  # nfeatures_val = 1500
  # data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val) 
  # 
  # pathName = paste0(folder_input,'QC/FindVariableFeatures',file_str,'.png')
  # print(pathName)
  # png(file=pathName,width=600, height=350, res = 100)
  # print(VariableFeaturePlot(data) + ylim(0,10))
  # dev.off()
  # 
  
  #Score for Cell Cycle gene expression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  #browser()
  # Scale
  # Do we need to scale?
  if (scale_TF){ # Check if we need to regress out
    # CHANGE vars.to.regress WITH VALUES FROM PARAM
    data <- ScaleData(data, features = rownames(data), vars.to.regress = c("nCount_RNA", "percent.mt"))
  }else{
    data <- ScaleData(data, features = rownames(data)) 
  }
  # PCA
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  data = RunPCA(data, features = VariableFeatures(object = data), npcs = PCA_dim)

  data = visualize_PCA(data,PCA_dim)
  #JackStrawPlot(data, dims = 1:PCA_dim)
  
  
  pathName <- paste0(folder_input,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  
  
  # Cluster with Umap
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  print(paste0('Resolution' ,': ', resolution_val))
  print(folder_input)
  data = getCluster (data,resolution_val, PCA_dim)
  
  # Name cells
  #data = label_cells(data,folder_input,sample_name,sampleParam,resolution_val,filter,scale_TF, file_str)
  
  ###############################
  # Find Cluster Biomarkers
  ###############################
  print('Find Cluster Biomarkers')
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

  
  all_markers =  markers %>% group_by(cluster)
  
  # Add known markers
  all_markers = cellMarkers(all_markers,cell_features)
  print(paste0(folder_input,'DE/Top20Features',file_str,'_Patient',Patient_num,'.csv'))
  write.csv(all_markers, file = paste0(folder_input,'DE/AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  
  # Get DE between cluster permutations
  filepath_mvn = paste0( folder_input, 'DE/nVsm/')
  dir.create( filepath_mvn, recursive = TRUE)
  Features_mvn_df = diffExpPermute(data,cell_features)
  
  level_list = levels(data)
  for (level in 1:length(levels(data)) ){
    ident1 = level_list[level]
    
    df_output = Features_mvn_df[[level]]
    df_output = df_output[c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","ident_1","ident_2","gene", "Cell")]
    write.csv(df_output, file = paste0(filepath_mvn
                                       ,'Features_',ident1,'Vsn',file_str
                                       ,'_Patient',Patient_num,'.csv'),row.names = FALSE)
  }  

  ########################
  ## Visualize clustering
  ########################
  filepath_cluster = paste0( folder_input, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  dir.create( filepath_cluster, recursive = TRUE)
  print(filepath_cluster)
  
  # Umap
  pathName = paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  # Cluster Metrics
  pathName <- paste0(filepath_cluster,paste0('ClusterMetrics',file_str,'.png'))
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  # Cluster Heatmap
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data, features = top10$gene))
  dev.off()
  
  return(data)
  
}