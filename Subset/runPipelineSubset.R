runPipelineSubset = function(data,folder_input,sample_name,sampleParam,filter,regress_TF){
  
  require(gtools)
  print(sample_name)
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
  
  cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
  
  file_str = ''
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)

  cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features <- read_excel(cell_features_file)
  
  # QC
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  print('hello')
  #data = quality_control(data,filter,nFeature_RNA_list,percent_mt,sample_name)
  
  #Score for Cell Cycle gene expression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  

  # Scale
  if (regress_TF){ # Check if we need to regress out
    # CHANGE vars.to.regress WITH VALUES FROM PARAM
    data <- ScaleData(data, features = rownames(data), vars.to.regress = c("nCount_RNA", "percent.mt"))
  }else{
    data <- ScaleData(data, features = rownames(data)) 
  }
  
  # PCA
  PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  data <- RunPCA(data, features = VariableFeatures(object = data), npcs = PCA_dim)
  visualize_PCA(data,folder_input,sample_name,PCA_dim)
  
  data = visualize_dim(data,PCA_dim)
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
  #data = label_cells(data,folder_input,sample_name,sampleParam,resolution_val,filter,regress_TF, file_str)
  
  
  
  
  ########################
  ## Visualize clustering
  ########################
  #browser()
  filepath_cluster = paste0( folder_input, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  dir.create( filepath_cluster, recursive = TRUE)
  print(filepath_cluster)
  # Plot Umap
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
  
  
  ###############################
  # Find Cluster Biomarkers
  ###############################
  print('Find Cluster Biomarkers')
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top20 = markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  
  
  all_markers =  markers %>% group_by(cluster)
  
  # Add known markers
  top20 = cellMarkers(top20,cell_features)
  all_markers = cellMarkers(all_markers,cell_features)
  print(paste0(folder_input,'DE/Top20Features',file_str,'_Patient',Patient_num,'.csv'))
  write.csv(top20, file = paste0(folder_input,'DE/Top20Features',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  write.csv(all_markers, file = paste0(folder_input,'DE/AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  
  # Get DE between cluster permutations
  perm = permutations(n = length(levels(data)), r = 2, v = 0:length(levels(data)))
  #perm = t(apply(perm,1,sort))
  #perm = unique( perm )
  
  filepath_mvn = paste0( folder_input, 'DE/nVsm/')
  dir.create( filepath_mvn, recursive = TRUE)
  
  Features_mvn_output = data.frame(matrix(ncol = 9, nrow = 0))
  names(Features_mvn_output) = c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","gene", "Cell","ident_1","ident_2")
  Features_mvn_df = vector(mode = "list", length = length(levels(data)))
  
  
  for (level in 1:length(levels(data)) ){
    Features_mvn_df[[level]] = Features_mvn_output
  }
  
  for (i in 1:nrow(perm)){
    
    
    Features_mvn = FindMarkers(data, ident.1 = perm[i,1], ident.2 = perm[i,2]
                               ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = TRUE)
    Features_mvn$gene = rownames(Features_mvn)
    rownames(Features_mvn) = NULL
    Features_mvn = cellMarkers(Features_mvn,cell_features)
    
    data_level = levels(data)
    Features_mvn$ident_1 = data_level[perm[i,1] + 1]
    Features_mvn$ident_2 = data_level[perm[i,2] + 1]
    
    Features_mvn_df[[perm[i,1] + 1]] = rbind(Features_mvn_df[[perm[i,1] + 1]], Features_mvn)
  }
  for (level in 1:length(levels(data)) ){
    df_output = Features_mvn_df[[level]]
    df_output = df_output[c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","ident_1","ident_2","gene", "Cell")]
    write.csv(df_output, file = paste0(filepath_mvn
                                       ,'Features_',(level-1),'Vsn',file_str
                                       ,'_Patient',Patient_num,'.csv'),row.names = FALSE)
  }
  
  return(data)
  
}