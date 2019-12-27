IntegrateAll_ClusterUmap = function(data,sample_type,folder_base_output,PCA_dim,resolution_val,label){
  #browser()
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
    
    plotHLA(data_dexaT,folder_base_output,paste0(sample_type,'_dexaT',label_str))
    plotHLA(data_dexaF,folder_base_output,paste0(sample_type,'_dexaF',label_str))
    
  }else if (sample_type == 'PreNBA'|| sample_type == 'PostNBM'){
    
    data_dexaT = SubsetData(object = data, cells = which(data$dexa == "Yes"))
    
    data_dexaF = SubsetData(object = data, cells = which(data$dexa == "No"))
    
    data_NBM = SubsetData(object = data, cells = which(data$dexa == "NBM"))
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_',sample_type,label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaT, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_',sample_type,label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_dexaF, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_NBM_',sample_type,label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_NBM, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    plotHLA(data_dexaT,folder_base_output,paste0(sample_type,'_dexaT',label_str))
    plotHLA(data_dexaF,folder_base_output,paste0(sample_type,'_dexaF',label_str))
    plotHLA(data_NBM,folder_base_output,paste0(sample_type,'NBA',label_str))
    
  }else if (sample_type == 'PrePost'){
    
    #data_post =SubsetData(object = data, cells = which())  
    
    data_dexaT_pre = SubsetData(object = data, cells = (data$split_var == "Pre D"))
    data_dexaF_pre = SubsetData(object = data, cells = (data$split_var == "Pre ND"))
    
    data_dexaT_post = SubsetData(object = data, cells =(data$split_var == "Post D"))
    data_dexaF_post = SubsetData(object = data, cells = (data$split_var == "Post ND"))
    
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
    
    plotHLA(data_dexaT_pre,folder_base_output,paste0(sample_type,'_dexaT_pre',label_str))
    plotHLA(data_dexaF_pre,folder_base_output,paste0(sample_type,'_dexaF_pre',label_str))
    plotHLA(data_dexaT_post,folder_base_output,paste0(sample_type,'_dexaT_post',label_str))
    plotHLA(data_dexaF_post,folder_base_output,paste0(sample_type,'_dexaF_post',label_str))
    
    
  }else if (sample_type == 'PrePostNBM'){
    
    data_dexaT_pre = SubsetData(object = data, cells = (data$split_var == "Pre D"))
    data_dexaF_pre = SubsetData(object = data, cells = (data$split_var == "Pre ND"))
    
    data_dexaT_post = SubsetData(object = data, cells =(data$split_var == "Post D"))
    data_dexaF_post = SubsetData(object = data, cells = (data$split_var == "Post ND"))
    
    data_NBM = SubsetData(object = data, cells = (data$dexa == "NBM"))
    
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
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_NBM_','Post',label_str,'.png'))  
    png(file=pathName,width=600, height=350)
    print(DimPlot(data_NBM, label=T, repel=F, reduction = "umap"))
    dev.off()
    
    
    plotHLA(data_dexaT_pre,folder_base_output,paste0(sample_type,'_dexaT_pre',label_str))
    plotHLA(data_dexaF_pre,folder_base_output,paste0(sample_type,'_dexaF_pre',label_str))
    plotHLA(data_dexaT_post,folder_base_output,paste0(sample_type,'_dexaT_post',label_str))
    plotHLA(data_dexaF_post,folder_base_output,paste0(sample_type,'_dexaF_post',label_str))
    
    plotHLA(data_NBM,folder_base_output,paste0(sample_type,'_NBM'))
  }
  
  
}