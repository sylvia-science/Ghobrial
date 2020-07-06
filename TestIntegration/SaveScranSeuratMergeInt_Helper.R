
SaveScranSeuratMergeInt_Helper = function(sample_list,mode,
                                          folder_output_merge,folder_output_integrate,
                                          downsample = NA){
  
  browser()
  Umap_type =  'PCA_Umap'
  dir.create(folder_output_merge,recursive = T)
  dir.create(folder_output_integrate,recursive = T)
  
  data_list = vector(mode = "list",length = 0)
  data_norm_list = data_list
  
  for (i in 1:length(sample_list)){
    #browser()
    sample_name = sample_list[i]
    print(sample_name)
    folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name , '/')
    data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
    data_i$sample = sample_name
    path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
    cellIdents = read.csv(path,sep = ',',row.names = 1)
    cellIdents$x = paste0(cellIdents$x, ' S',i)
    data_i$CellType = cellIdents
    
    if (!is.na(downsample)){
      downsample = sub("_.*", "", downsample)
      cellnames = colnames(data_i)
      cellnames = sub("_.*", "", cellnames)
      data_i = data_i[,cellnames %in% downsample]
      #browser()
    }

    print(ncol(data_i))
    if (ncol(data_i) > 100){
      data_list = c(data_list, data_i)
        
      data_i_norm = ScranNorm(data_i)
      data_i_norm <- FindVariableFeatures(data_i_norm, selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
        
      data_norm_list = c(data_norm_list,data_i_norm)
    }
    
  
  }
  
  
  browser()
  data_norm_list =data_norm_list[lengths(data_norm_list) != 0]
  if (mode == 'Merge'){
    
    data_merge = merge(x = data_norm_list[[1]] ,y = data_norm_list[2:length(data_norm_list)], merge.data = T)
    #for (i in 3:length(data_norm_list)){
    #  print(i)
    #  data_merge = merge(x =  data_merge,y = data_norm_list[[i]],merge.data = T)
    #}
    
    #browser()
    resolution_val = 1.4
  
  
    
    #data_merge_run = ScranNorm(data_merge)
    data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
    data_merge_run = ScaleData(data_merge_run)
    data_merge_run = RunPCA(data_merge_run,npcs = 30)
    data_merge_run = FindNeighbors(data_merge_run, dims = 1:30)
    data_merge_run = FindClusters(data_merge_run, resolution = resolution_val)
    #data_merge_run = RunUMAP(data_merge_run,umap.method = 'umap-learn',graph= 'RNA_snn')
    
    if (Umap_type == 'PCA_Umap'){
      data_merge_run = RunUMAP(data_merge_run, dims = 1:30)
    }else{
      data_merge_run = RunUMAP(data_merge_run,umap.method = 'umap-learn',graph= 'RNA_snn')
    }
    
    data_merge_run$sample = data_merge$sample
    png(file=paste0(folder_output_merge,'merge.png'),width=1000, height=1000)
    plot <- DimPlot(data_merge_run, reduction = "umap", group.by = "sample")
    print(plot)
    dev.off()
    
    data_merge_run$split_var = ''
    cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
    groupBy_list = c('sample')
    splitBy_list = c('sample')
    plotAll(data_merge_run, folder = folder_output_merge,
            sample_name,sampleParam = NA,
            cell_features = cell_features,
            label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
            clusterTF =F, markersTF = F, keepOldLabels = F, 
            groupBy = groupBy_list, splitBy = splitBy_list,
            PCA_dim = 30,resolution_val = resolution_val)
    
    filepath_cluster = paste0( folder_output_merge, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
    #PlotKnownMarkers(data_merge_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
    #                 plotType ='FeaturePlotFix' , str = '')
  
   
    #browser()
    
    path = paste0(filepath_cluster,'','data','.Robj')
    save(data_merge_run,file= path)
    #next
  }else if (mode == 'Integrate'){
    ######################################
    ## Integrate
    
   
    browser()
    k.anchor = 5
    k.filter = 100
    k.score = 30
    k.weight = 200
    anchors = FindIntegrationAnchors(object.list = data_norm_list,
                                     k.anchor = k.anchor, k.filter =k.filter, k.score = k.score, 
                                     anchor.features = 2000)
    #browser()
    resolution_val = 1.4
    filepath_cluster = paste0( folder_output_integrate, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
    dir.create(filepath_cluster, recursive = T)
    path = paste0(filepath_cluster,'','anchors',
                  'kanchor_',k.anchor,'kfilter',k.filter,'kscore',k.score,'.Robj')
    save(anchors,file= path)
    anchors = loadRData(path)

    data_integrate = IntegrateData(anchorset = anchors, 
                                   dims = 1:30,k.weight  = k.weight)
    #browser()
    
    path = paste0(filepath_cluster,'','data_integrate','.Robj')
    save(data_integrate,file= path)
    
    path = paste0(filepath_cluster,'','data_integrate','.Robj')
    data_integrate  = loadRData(path)
    
    data_integrate <- ScaleData(data_integrate, verbose = FALSE)
    data_integrate <- RunPCA(data_integrate, npcs = 30, verbose = FALSE)
    data_integrate = FindNeighbors(data_integrate, dims = 1:30)
    data_integrate = FindClusters(data_integrate, resolution = resolution_val)
    #browser()
    
    if (Umap_type == 'PCA_Umap'){
      data_integrate = RunUMAP(data_integrate, dims = 1:30)
    }else{
      data_integrate <-  RunUMAP(data_integrate,umap.method = 'umap-learn',graph= 'integrated_snn')
    }
    
    path = paste0(filepath_cluster,'','data','.Robj')
    save(data_integrate,file= path)
    
    #################
    
    data_integrate$split_var = ''
    cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
    groupBy_list = c('sample')
    splitBy_list = c('sample')
    plotAll(data_integrate, folder = folder_output_integrate,
            sample_name,sampleParam = NA,
            cell_features = cell_features,
            label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
            clusterTF =F, markersTF = F, keepOldLabels = F, 
            groupBy = groupBy_list, splitBy = splitBy_list,
            PCA_dim = 30,resolution_val = resolution_val)
    
    filepath_cluster = paste0( folder_output_integrate, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
    PlotKnownMarkers(data_integrate, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                     plotType ='FeaturePlotFix' , str = '')
    
    
  
  }
}