run_pipeline = function(data,folder_input,sample_name,sampleParam,filter,regress_TF){
  
  #browser()
  remove_NaN = F
  if (remove_NaN){
    # Remove cells with Nan values
    tmp = data@assays[["RNA"]]@data
    tmp = as.matrix(tmp)
    cell_name = names(which(is.na(colSums(tmp))))
    cell_val = colnames(data[["RNA"]]@data) != cell_name
    colname_list = colnames(data[["RNA"]]@data)
    # Stupid fix
    for ( i in 1:length(colname_list)){
      if(colname_list[i] %in% cell_name) {
        cell_val[i] = F
      }
      
    }
    
    #data_integrate_test = SubsetData(data, cells = cell_val)
    if (length(cell_val)!= 0){
      data = data[,cell_val]
    }
  }
  #browser()
  dir.create(folder_input,recursive = TRUE)
  print(sample_name)
  #browser()
  cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  
  file_str = ''
  
  
  cell_features_file = '/home/sujwary/Desktop/scRNA/Data/Cell_IDS.xlsx'
  cell_features = read_excel(cell_features_file)
  #browser()
  #####################
  ## QC
  #####################
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  folder_input = makeFolders(folder_input,sample_name,filter,regress_TF, makeFolder_TF= TRUE,nFeature_RNA_list,percent_mt)
  
  print('hello')
  data = quality_control(data,folder_input,filter,nFeature_RNA_list,percent_mt,sample_name)
  data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  ########################
  # Get Variable Genes
  ########################
  
  nfeatures_val = sampleParam$nfeatures_val[sampleParam['Sample'] == sample_name]
  data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val)
  
  pathName = paste0(folder_input,'QC Metrics/FindVariableFeatures',file_str,'.png')
  print(pathName)
  png(file=pathName,width=600, height=350, res = 100)
  print(VariableFeaturePlot(data) + ylim(0,10))
  dev.off()
  
  ##########################################
  ## Score for Cell Cycle gene expression
  ##########################################
  s.genes = cc.genes$s.genes
  g2m.genes = cc.genes$g2m.genes
  data = CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
  # Scale
  if (regress_TF){ # Check if we need to regress out
    # CHANGE vars.to.regress WITH VALUES FROM PARAM
    data <- ScaleData(data, features = rownames(data), vars.to.regress = c("nCount_RNA", "percent.mt"))
  }else{
    data <- ScaleData(data, features = rownames(data)) 
  }
  #browser()
  # PCA
  PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  data <- RunPCA(data, features = VariableFeatures(object = data), npcs = PCA_dim)
  visualize_PCA(data,folder_input,PCA_dim)
  
  #data = visualize_dim(data,PCA_dim)
  #JackStrawPlot(data, dims = 1:PCA_dim)
  
  
  pathName <- paste0(folder_input,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  
  
  # Cluster with Umap
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  print(paste0('Resolution' ,': ', resolution_val))
  print(folder_input)
  data <- getCluster (data,resolution_val, PCA_dim)
  
  # Name cells
  #data = label_cells(data,folder_input,sample_name,sampleParam,resolution_val,filter,regress_TF, file_str)
  
  
  
  # Find Cluster Biomarkers
  print('Find Cluster Biomarkers')
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  
  
  
  filepath_cluster = paste0( folder_input, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  # 
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data, features = top10$gene,label = F))
  dev.off()
  
  #browser()
  
  
  # Add known markers
  all_markers =  markers %>% group_by(cluster)
  all_markers = cellMarkers(all_markers,cell_features)
  write.csv(all_markers, file = paste0(filepath_cluster,'Features',file_str,'.csv'),row.names=FALSE)
  #browser()
  
  
  #browser()
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data, features = top10$gene, label=FALSE))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  #########################################################################################
  
  
  return(data)
  
}