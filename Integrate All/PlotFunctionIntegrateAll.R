IntegrateAll_ClusterUmap = function(data,sample_type,folder_base_output,PCA_dim,resolution_val,label){
  
  filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  if (label){
    label_str = '_label'
  }else{
    label_str = ''
  }
  
  if (sample_type == 'Pre' || sample_type == 'Post'){
    
    
    # Grab dexa == 'Yes'
    data_dexaT = SubsetData(object = data, cells = which(data$dexa == "Yes"))
    
    # Grab dexa == 'No'
    data_dexaF = SubsetData(object = data, cells = which(data$dexa == "No"))
    
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_',sample_type,label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaT, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_',sample_type,label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaF, label=T, repel=F, reduction = "umap"))
    dev.off()
    
  }else if (sample_type == 'PrePost'){
    
    data_post =SubsetData(object = data, cells = which())  
    
    # Grab dexa == 'Yes'
    data_dexaT_pre = SubsetData(object = data, cells = which(data$dexa == "Yes" && data$orig.ident == "data_pre"))
    
    # Grab dexa == 'No'
    data_dexaF_pre = SubsetData(object = data, cells = which(data$dexa == "No"&& data$orig.ident == "data_pre"))
    
    
    # Grab dexa == 'Yes'
    data_dexaT_post = SubsetData(object = data, cells = which(data$dexa == "Yes" && data$orig.ident == "data_post"))
    
    # Grab dexa == 'No'
    data_dexaF_post = SubsetData(object = data, cells = which(data$dexa == "No"&& data$orig.ident == "data_post"))
    
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_','Pre',label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaT_pre, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_','Pre',label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaF_pre, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_','Post',label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaT_post, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_','Post',label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaF_post, label=T, repel=F, reduction = "umap"))
    dev.off()
  }
  
}

#########################################
