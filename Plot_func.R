plotAll = function(data,folder,sample_name,sampleParam,label_TF){
  print(sample_name)
  print(folder)
  print(paste0('label_TF: ', label_TF))
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]

  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)

  cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features <- read_excel(cell_features_file)

  #cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  if (label_TF){
    file_str = '_label'
  }else
  {
    file_str = ''
  }
  print(paste0('file_str: ', file_str))

  # PCA
  PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  visualize_PCA(data,folder,sample_name,PCA_dim)
  visualize_dim(data)

  pathName <- paste0(folder,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  # Cluster with Umap
  ## ADD DOTPLOT
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  data <- getCluster (data,resolution_val, PCA_dim)
  
  
  # Name cells
  if (label_TF){
    
    cluster_IDs <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
    #browser()
    data = label_cells(data,cluster_IDs)
  }
  
  
  ##############################
  ## Find Cluster Biomarkers
  ##############################
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  all_markers =  markers %>% group_by(cluster)
  
  # Add known markers to top20
  top20 = cellMarkers(top20,cell_features)
  all_markers = cellMarkers(all_markers,cell_features)
  
  write.csv(top20, file = paste0(folder,'DE/Top20Features',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  write.csv(all_markers, file = paste0(folder,'DE/AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  
  ########################
  ## Visualize clustering
  ########################
  pathName <- paste0(folder,paste0('Cluster/ClusterUmap',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(folder,paste0('Cluster/HeatMap',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data, features = top10$gene))
  dev.off()
  
  pathName <- paste0(folder,'Cluster/ClusterMetrics','.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()

  
}

##
## Integrate
##
# 
# plotAllIntegrate = function(data,folder,sample_name,sampleParam,label_TF){
#   print(sample_name)
#   print(folder)
#   print(paste0('label_TF: ', label_TF))
#   
#   Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
# 
#   # Load data
#   filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
#   metaData <- read_excel(filename_metaData)
# 
#   cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
#   cell_features <- read_excel(cell_features_file)
#   
#   #cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
#   if (label_TF){
#     file_str = '_label'
#   }else
#   {
#     file_str = ''
#   }
#   
#   print(paste0('file_str: ', file_str))
#   
#   # PCA
#   PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
#   visualize_PCA(data,folder,sample_name,PCA_dim)
#   visualize_dim(data)
#   
#   pathName <- paste0(folder,'PCA/elbow_',PCA_dim,'.png')
#   png(file=pathName,width=600, height=350)
#   print(ElbowPlot(data,ndims = PCA_dim))
#   dev.off()
#   
#   
#   
#   # Cluster with Umap
#   ## ADD DOTPLOT
#   resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
#   data <- getCluster (data,resolution_val, PCA_dim)
#   
#   
#   # Name cells
#   if (label_TF){
#     cluster_IDs <- sampleParam$Cluster_IDs_BM[sampleParam['Sample'] == sample_name]
#     data = label_cells(data,cluster_IDs)
#   }
#   
#   ##############################
#   ## Find Cluster Biomarkers
#   ##############################
#   # find markers for every cluster compared to all remaining cells, report only the positive ones
#   markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#   markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#   # Plotting the top 10 markers for each cluster.
#   top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#   top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#   all_markers <- markers %>% group_by(cluster)
#   
#   top20$Cell = NA
#   # Add known markers 
#   
#   top20 = cellMarkers(top20,cell_features)
#   all_markers = cellMarkers(all_markers,cell_features)
#   write.csv(top20, file = paste0(folder,'Top20Features',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
#   write.csv(all_markers, file = paste0(folder,'AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
#   
#   ########################
#   ## Visualize clustering
#   ########################
#   pathName <- paste0(folder,paste0('Cluster/ClusterUmap',resolution_val,file_str,'.png'))
#   png(file=pathName,width=600, height=350, res = 100)
#   print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
#   dev.off()
#   
#   pathName <- paste0(folder,paste0('Cluster/HeatMap',resolution_val,file_str,'.png'))
#   png(file=pathName,width=1000, height=1200)
#   print(DoHeatmap(data, features = top10$gene))
#   dev.off()
#   
#   pathName <- paste0(folder,'Cluster/ClusterMetrics','.png')
#   png(file=pathName,width=600, height=350)
#   print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
#   dev.off()
# 
#   # Split by orig ident
#   pathName <- paste0(folder,paste0('Cluster/ClusterUmap',resolution_val,'_split',file_str,'.png'))  
#   png(file=pathName,width=600, height=350)
#   print(DimPlot(data, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
#   dev.off()
#   
#   return(data)
#   
# }
