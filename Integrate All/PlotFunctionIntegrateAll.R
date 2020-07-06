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
FeaturePlot_GeneList = function(data,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = FALSE){
  #browser()
  library(cowplot)
  folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
  dir.create( folder_featureplot, recursive = TRUE)
  print(folder_featureplot)
  for (i in 1:(length(gene_list))){
    feature = as.character(gene_list[i])
    print(feature)
    #FeaturePlotFix(data, feature,folder_featureplot,'', split = FALSE, gene_TF = TRUE,title = folder_name)
    #MakeFeaturePlot(data,data, folder_featureplot, feature,split = FALSE)
    
    if (!FeaturePlotFix){
    
      pathName <- paste0(folder_featureplot,paste0( feature,'.png'))  
      png(file=pathName,width=1500, height=600)
      print(FeaturePlot(data, features = feature,split.by = "split_var"))
      dev.off()
    }else{

      label_str = ''
      if (sample_type == 'Pre' || sample_type == 'Post'){
        
        data_dexaT = SubsetData(object = data, cells = (data$dexa == "Yes"))
        data_dexaF = SubsetData(object = data, cells = (data$dexa == "No"))
        
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_',sample_type,label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaT, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_',sample_type,label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaF, label=T, repel=F, reduction = "umap"))
        dev.off()
        
      }else if (sample_type == 'PreNBA'|| sample_type == 'PostNBM'){
        
        data_dexaT = SubsetData(object = data, cells = (data$dexa == "Yes"))
        
        data_dexaF = SubsetData(object = data, cells = (data$dexa == "No"))
        
        data_NBM = SubsetData(object = data, cells = (data$dexa == "NBM"))
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_',sample_type,label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaT, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_',sample_type,label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaF, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_NBM_',sample_type,label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_NBM, label=T, repel=F, reduction = "umap"))
        dev.off()

      }else if (sample_type == 'PrePost'){
        
        #data_post =SubsetData(object = data, cells = which())  
        
        data_dexaT_pre = SubsetData(object = data, cells = (data$split_var == "Pre D"))
        data_dexaF_pre = SubsetData(object = data, cells = (data$split_var == "Pre ND"))
        
        data_dexaT_post = SubsetData(object = data, cells =(data$split_var == "Post D"))
        data_dexaF_post = SubsetData(object = data, cells = (data$split_var == "Post ND"))
        
        #browser()
        plot_PreD = FeaturePlotFix(data_dexaT_pre, feature,folder = '','', split = FALSE, gene_TF = TRUE,title = '',saveTF = FALSE)
        plot_PreND = FeaturePlotFix(data_dexaF_pre, feature,folder = '','', split = FALSE, gene_TF = TRUE,title = '',saveTF = FALSE)
        plot_PostD = FeaturePlotFix(data_dexaT_post, feature,folder = '','', split = FALSE, gene_TF = TRUE,title = '',saveTF = FALSE)
        plot_PostND = FeaturePlotFix(data_dexaF_post, feature,folder = '','', split = FALSE, gene_TF = TRUE,title = '',saveTF = FALSE)
        
        label_list = c('Pre ND','Pre D','Post ND','Post D')
        plot = plot_grid(plot_PreND, plot_PreD,plot_PostND,plot_PostD,labels=label_list)

        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_',label_str,feature,'.png'))  
        png(file=pathName,width=1200, height=600)
        print(plot)
        dev.off()
        
      }else if (sample_type == 'PrePostNBM'){
        
        data_dexaT_pre = SubsetData(object = data, cells = (data$split_var == "Pre D"))
        data_dexaF_pre = SubsetData(object = data, cells = (data$split_var == "Pre ND"))
        
        data_dexaT_post = SubsetData(object = data, cells =(data$split_var == "Post D"))
        data_dexaF_post = SubsetData(object = data, cells = (data$split_var == "Post ND"))
        
        data_NBM = SubsetData(object = data, cells = (data$dexa == "NBM"))
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_','Pre',label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaT_pre, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_','Pre',label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaF_pre, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaT_','Post',label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaT_post, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_DexaF_','Post',label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_dexaF_post, label=T, repel=F, reduction = "umap"))
        dev.off()
        
        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_NBM_','Post',label_str,'.png'))  
        png(file=pathName,width=600, height=350)
        print(DimPlot(data_NBM, label=T, repel=F, reduction = "umap"))
        dev.off()
        
      
      }else if (sample_type == ''){
        plot = FeaturePlotFix(data, feature,folder = '','', split = FALSE, gene_TF = TRUE,title = '',saveTF = FALSE)

        pathName <- paste0(folder_featureplot,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_',label_str,feature,'.png'))  
        png(file=pathName,width=1200, height=600)
        print(plot)
        dev.off()
        
      }
    }
  }
}