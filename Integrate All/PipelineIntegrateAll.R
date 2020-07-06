PipelineIntegrateAll = function(data,
                                sample_name,
                                folder_output,
                                sampleParam,
                                integrate_merge,
                                regress_var,
                                markersTF = TRUE,
                                cell_features,
                                ConvertCategorical){
  file_str = ''
  
  
  browser()
  # Remove cells with Nan values
  tmp = data[[data_integrate@active.assay]]@data
  tmp = as.matrix(tmp)
  cell_name = names(which(is.na(colSums(tmp))))
  cell_val = colnames(data[[data_integrate@active.assay]]@data) != cell_name
  colname_list = colnames(data[[data_integrate@active.assay]]@data)
  # Stupid fix
  for ( i in 1:length(colname_list)){
    if(colname_list[i] %in% cell_name) {
      cell_val[i] = F
    }
    
  }
  #data_integrate_test = SubsetData(data, cells = cell_val)
  if (length(cell_val) !=0){
    data = data[,cell_val]
  }
  #browser()
  
  
  library(stringr)
  
  if (sample_name == 'Integrate_PrePostNBM_filterF' ){
    data_new = CreateSeuratObject(
      data@assays$RNA@data,
      assay = "RNA",
      min.cells = 1,
    )
    
    data_new@meta.data = data@meta.data
    data = data_new
  }
  #data@assays$RNA@data = data_new@assays$RNA@data
  # Load data
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample' ] == sample_name]
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  
  if (is.na(PCA_dim) || is.na(resolution_val) ){
    print('No PCA/Resolution found!')
    browser()
  }
  filepath_cluster = paste0( folder_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
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
    data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    
    ########################
    # Get Variable Genes
    ########################
    
    nfeatures_val =  2000
    data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val)
    
    dir.create( paste0(folder_output,'QC Metrics/'), recursive = TRUE)
    
    pathName = paste0(folder_output,'QC Metrics/FindVariableFeatures','.png')
    print(pathName)
    png(file=pathName,width=600, height=350, res = 100)
    print(VariableFeaturePlot(data) + ylim(0,10))
    dev.off()
    
    
  }
  
  #nfeatures_val = 2000 
  #data = FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val)

  
  #Score for cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  #browser()
  print(paste0('PCA: ', PCA_dim))
  
  #data_new <- CreateSeuratObject(data@assays$RNA@counts)
  #data_new$orig.ident <- data$orig.ident
  #data=  data_new 
  data = CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
  #browser()
  #data_prescale = data
  #regress_var = c(regress_var,c("nCount_RNA", "percent.mt") )
  if (regress_var!= ''){    
    if (ConvertCategorical == 'ConvCatT'){
      ## Convert variables to numeric for regression
      for (var in regress_var){
        data@meta.data[,var] = as.character(factor(data@meta.data[,var]))
        data_var  = unique(data@meta.data[,var])
        cnt = 1
        data@meta.data[,paste0(var,'_numeric')] =  NA
        for (level in data_var){
          data@meta.data[,paste0(var,'_numeric')][data@meta.data[,var] == level] = as.numeric(cnt)
          cnt = cnt + 1
        }
        
      }
      #browser()
      
      if ('PatientLast' %in% regress_var){
        browser()
        regress_var = regress_var[!str_detect(regress_var,pattern="PatientLast")]
        data = ScaleData(data, vars.to.regress = paste0(regress_var,'_numeric'))
        data = ScaleData(data, vars.to.regress = paste0('Patient','_numeric'))
      }else{
        data = ScaleData(data, vars.to.regress = paste0(regress_var,'_numeric'))
      }
    }else{
      #browser()
      if ('PatientLast' %in% regress_var){
        regress_var = regress_var[!str_detect(regress_var,pattern="PatientLast")]
        data = ScaleData(data, vars.to.regress = paste0(regress_var))
        data = ScaleData(data, vars.to.regress = paste0('Patient'))
      }else{
        #data = ScaleData(data)
        
        data = ScaleData(data, vars.to.regress = paste0(regress_var))
      }
    }
    
    
  }else{
    data = ScaleData(data)
    #data = ScaleData(data,vars.to.regress = c("nCount_RNA", "percent.mt"))
  }
  
  
  #data = FindVariableFeatures(data_integrate_test, selection.method = "vst", nfeatures = 2000)
  #data = data_integrate_test
  #browser()
  if (integrate_merge == 'Integrate'){
    
    #data_prePCA = data
    data = RunPCA(data, npcs = PCA_dim)
  }else{
    #data_prePCA = data
    data = RunPCA(data, features = data@assays[["RNA"]]@var.features, npcs = PCA_dim)
  }
  #browser()
  ######################
  ## Start Plotting
  ######################
  
  #Visualize PCA results
  visualize_PCA(data,folder_output,PCA_dim)
  
  pathName = paste0(folder_output,'PCA/elbow',PCA_dim,file_str,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  # Find Neighbors and Clusters
  data = FindNeighbors(data, dims = 1:PCA_dim)
  data = FindClusters(data, resolution = resolution_val) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.)
  data = RunUMAP(data, dims = 1:PCA_dim)
  
  # Visualize clustering
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1000, res = 100)
  print(DimPlot(data,pt.size = 0.5, reduction = "umap",label = TRUE))
  dev.off()
  
  # pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  # png(file=pathName,width=600, height=350)
  # print(FeaturePlot(data,pt.size = 0.5, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  # dev.off()
  
  if (markersTF == TRUE){
    # Find Cluster Biomarkers
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    # Plotting the top 10 markers for each cluster.
    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    all_markers =  markers %>% group_by(cluster)
    
    
    pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
    png(file=pathName,width=1000, height=1200)
    print(DoHeatmap(data, features = top10$gene))
    dev.off()
    
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
    
    
  }
  
  
  
  
  
  return(data)
  
}