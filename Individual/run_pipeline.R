run_pipeline = function(filename,folder_input,sample_name,sampleParam,filter,regress_TF){
  print(sample_name)
  
  cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  
  file_str = ''
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  metaData <- read_excel(filename_metaData)
  sampleParam <- read_excel(filename_sampleParam)
  data = load_data(filename)
  
  cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features <- read_excel(cell_features_file)
  
  #####################
  ## QC
  #####################
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  print('hello')
  data = quality_control(data,filter,nFeature_RNA_list,percent_mt,folder_input,sample_name)
  
  ########################
  # Get Variable Genes
  ########################
  nfeatures_val = sampleParam$nfeatures_val[sampleParam['Sample'] == sample_name]
  data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val)

  pathName = paste0(folder_input,'QC/FindVariableFeatures',file_str,'.png')
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
  
  top20$Cell = NA
  # Add known markers to top20
  #browser()
  
  for (i in 1:nrow(top20)){
    marker_rows = grep(top20$gene[i], cell_features$Markers, value=TRUE)
    marker_idx = which( cell_features$Markers %in% marker_rows)
    #browser()
    if (length(marker_idx) > 0){
      #browser()
      cell_list = cell_features$Cell[marker_idx]
      cell_list = paste(cell_list, sep="", collapse=", ") 
      top20$Cell[i] = cell_list
    }else if (length(marker_idx) == 0 ){
      top20$Cell[i] = ''
    }
  }
  
  
  
  write.csv(top20, file = paste0(folder_input,'Top20Features',file_str,'.csv'),row.names=FALSE)
  #browser()
  # Plot Umap
  pathName <- paste0(folder_input,paste0('Cluster/ClusterUmap',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  # Visualize clustering
  # Cluster Metrics
  pathName <- paste0(folder_input,'Cluster/ClusterMetrics.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  # Cluster Heatmap
  plot = DoHeatmap(data, features = top10$gene)
  pathName <- paste0(folder_input,paste0('Cluster/HeatMap',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(plot)
  dev.off()
  
  
  # Get gene Descriptions
  #gene_desc_top10 =  get_gene_desc(top10)
  #gene_desc_top10
  #########################################################################################
  
  
  return(data)
  
}