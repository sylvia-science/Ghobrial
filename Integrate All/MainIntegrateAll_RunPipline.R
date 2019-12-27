# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/MainIntegrateAll_RunPipline.R')
rm(list = ls())
gc()

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/PipelineIntegrateAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/IntegrateAll_ClusterUmap.R')


library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

require(gridExtra)
require(data.table) 
integrate_merge = 'Merge'


sample_type_list = c('PrePostNBM','Pre','Post','PrePost','PreNBM','PostNBM','NBM')

sample_type_list = c('PreNBM','PostNBM','NBM')

sample_type_list = c('PrePost')

print(paste0('integrate_merge:', integrate_merge))

if (integrate_merge == 'Integrate' || integrate_merge == 'Merge'){
  filename_sampleParam <- paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_','Combine','_parameters.xlsx')
  sampleParam <- read_excel(filename_sampleParam)
}else{
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  sampleParam <- read_excel(filename_sampleParam)
}

print(filename_sampleParam)
i = 1
for (i in 1:length(sample_type_list)){
  sample_type = sample_type_list[i]
  sample_name = paste0(integrate_merge, '_',sample_type)
  folder_base_output = paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/',integrate_merge ,' All/',sample_type,'/')
  
  
  pathName <- paste0(folder_base_output,'PCA')
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder_base_output,'Cell Type')
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder_base_output,'Cluster')
  dir.create( pathName, recursive = TRUE)
  
  
  pathName <- paste0(folder_base_output,'DE')
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder_base_output,'Stats')
  dir.create( pathName, recursive = TRUE)
  

  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  
  
  RunPipeline = FALSE
  if (RunPipeline){
    print('Run')
    path = paste0(folder_base_output,'data','_',sample_name,'.Robj')
    print(path)
    data_integrate = loadRData(path)
    data_run = PipelineIntegrateAll(data_integrate,sample_name,folder_base_output,sampleParam,integrate_merge)
    #browser()  
    
    
    save(data_run,file=paste0(folder_base_output,'data_run_',integrate_merge,'_PCAdim',PCA_dim,'_',sample_type,'.Robj'))
    
    PlotKnownMarkers(data_run,data_run, paste0(folder_base_output,'Cell Type/'), cell_features = NA,split = FALSE)
    
    
  }else{
    #browser()
    print('Plot')
    path = paste0(folder_base_output,'data_run_',integrate_merge,'_PCAdim',PCA_dim,'_',sample_type,'.Robj')
    print(path)
    data_run = loadRData(path)
    data_run@meta.data$split_var = paste(data_run@meta.data$orig.ident,data_run@meta.data$dexa)
    data_run@meta.data$split_var = gsub("data_post Yes", "Post D", data_run@meta.data$split_var)
    data_run@meta.data$split_var = gsub("data_pre Yes", "Pre D", data_run@meta.data$split_var)
    data_run@meta.data$split_var = gsub("data_post No", "Post ND", data_run@meta.data$split_var)
    data_run@meta.data$split_var = gsub("data_pre No", "Pre ND", data_run@meta.data$split_var)
    data_run@meta.data$split_var = gsub("data_NBM NBM", "NBM", data_run@meta.data$split_var)
    
    
    
    #plotAll(data_run,folder_base_output,sample_name,sampleParam,label_TF = FALSE,integrate_TF = FALSE,  DE_perm_TF = FALSE,clusterTF = TRUE)
    #plotAll(data_run,folder_base_output,sample_name,sampleParam,label_TF = TRUE,integrate_TF = FALSE,  DE_perm_TF = FALSE,clusterTF = TRUE)
    #browser()
    #PlotKnownMarkers(data_run,data_run, paste0(folder_base_output,'Cell Type/'), cell_features = NA,split = FALSE)
    
    
    #browser()
    
    data_run = getCluster (data_run,resolution_val, PCA_dim)
    cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
    data_run_label = label_cells(data_run,cluster_IDs)
    
    
    gene_list = c('FCGR3A', 'SLAMF7', 'KLRK1', 'ICAM1', 'ITGAL', 'ITGB2','NCR1', 'SH2D1B','KLRC2', 'KLRC1', 'KLRD1')
    gene_list = c('GNLY','NKG7')
    gene_list = c('DUSP1','RRM2','BCL2L1','IKZF1','NR3C1')
    folder_name = 'Feature Plots 12-20-2019'
    gene_list = c('BIRC3', 'CFB', 'CSF3', 'CXCL2', 'EFNA1', 'G0S2', 'CSF2', 'IFIT1', 'IL32', 'LAMB3', 'PRIC285', 'TNFAIP3')
    
    #FeaturePlot_GeneList(data_run_label,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = FALSE)
    
    FeaturePlot_GeneList(data_run_label,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    
    browser()
    
    data_run_label = label_cells(data_run,gsub("dCD14+ Mono", "CD14+ Mono", cluster_IDs))
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
    
    ident1 = paste0('T Cell', ' ', 'Pre ND')
    ident2 = paste0('T Cell', ' ', 'Post ND')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)

    Features = Features[order(Features$avg_logFC),]
    
    gene_list = rownames(Features)
    gene_list = gene_list[1:50]
    folder_name = 'Dexa Gene Feature Plots 12-17-2019'
    FeaturePlot_GeneList(data_run_label,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    ##############
    ident1 = paste0('CD14+ Mono', ' ', 'Pre ND')
    ident2 = paste0('CD14+ Mono', ' ', 'Post ND')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    
    Features = Features[order(Features$avg_logFC),]
    Features = Features[Features$p_val_adj < 0.05,]
    gene_list = rownames(Features)
    gene_list = gene_list[1:50]
    folder_name = 'Dexa CD14+ ND Gene Feature Plots 12-17-2019'
    FeaturePlot_GeneList(data_run,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    
    #######################
    cluster_IDs_new = gsub("NK1", "NK", cluster_IDs)
    cluster_IDs_new = gsub("NK2", "NK", cluster_IDs_new)
    
    data_run_label = label_cells(data_run,cluster_IDs_new)
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
    
    ident1 = paste0('NK', ' ', 'Pre ND')
    ident2 = paste0('NK', ' ', 'Post ND')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    
    Features = Features[order(Features$avg_logFC),]
    Features = Features[Features$p_val_adj < 0.05,]
    gene_list = rownames(Features)
    gene_list = gene_list[1:50]
    folder_name = 'Dexa NK ND Gene Feature Plots 12-17-2019'
    FeaturePlot_GeneList(data_run,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    #######################
    
    
    gene_list = c('GZMA','GZMB','GZMH','IL2RB','TGFB1')
    folder_name = 'Dexa DNK ND Gene Feature Plots 12-17-2019'
    FeaturePlot_GeneList(data_run,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    #######################
    
    
    
    path = paste0(folder_base_output,'DE/','T','/')
    dir.create( path, recursive = TRUE)
    path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
    print(path)
    write.csv(Features, file = path,row.names=TRUE)
    
    
    browser()
    #PlotKnownMarkers(data_run,data_run, paste0(folder_base_output,'Cell Type/'), cell_features = NA,split = FALSE)
    # Label:
    
    #IntegrateAll_ClusterUmap(data_run,sample_type,folder_base_output,PCA_dim,resolution_val,label = FALSE)
    #IntegrateAll_ClusterUmap(data_run_label,sample_type,folder_base_output,PCA_dim,resolution_val,label = TRUE)
    
    
    #########################
    filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
    
    pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_splitAll', '','.png'))  
    png(file=pathName,width=2600, height=500,res = 100)
    print(DimPlot(data_run_label, label=T, repel=F, reduction = "umap", split.by = "split_var"))
    dev.off()
    ###############################
    
    #browser()
    

    
    #data_run = label_cells(data_run,cluster_IDs)
    
    
    # 
    # data = SubsetData(object = data_run, cells = data_run$split_var == 'Post D' )
    # plot = DimPlot(data, label=T, repel=F, reduction = "umap",pt.size = 1) + 
    #   ggtitle('Post D' ) + 
    #   theme(plot.title = element_text(hjust = 0.5))
    # pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_splitAll_PostD', '','.png'))  
    # png(file=pathName,width=1000, height=500,res = 100)
    # print(plot)
    # dev.off()
    # 
    # data = SubsetData(object = data_run, cells = data_run$split_var == 'Pre D' )
    # plot = DimPlot(data, label=T, repel=F, reduction = "umap", pt.size = 1) + 
    #   ggtitle('Pre D' ) + 
    #   theme(plot.title = element_text(hjust = 0.5))
    # pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_splitAll_PreD', '','.png'))  
    # png(file=pathName,width=1000, height=500,res = 100)
    # print(plot)
    # dev.off()

    
    #browser()
    
    ####################

    print('Get Markers')
    data_run_label = label_cells(data_run,gsub("dCD14+ Mono", "CD14+ Mono", cluster_IDs))
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
    
    cell_list = c('T Cell','Mono1','Mono2','Mono3', 'NK')
    category_list = c('Pre ND','Pre D', 'Post ND', 'Post D')
    
    for (cat_i in category_list){
      for (cat_j in category_list){
        for (cell in cell_list){
          if (cat_i != cat_j){
            ident1 = paste0(cell, ' ', cat_i)
            ident2 = paste0(cell, ' ', cat_j)
            Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                                                ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
            
            path = paste0(folder_base_output,'DE/',cell,'/')
            dir.create( path, recursive = TRUE)
            path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
            print(path)
            write.csv(Features, file = path,row.names=TRUE)
          }
        }
      }
      
    }

    

    ######################################33
    data = SubsetData(object = data_run_label, cells = (data_run$dexa == "Yes" && data_run$orig.ident == "data_post"))
    data_run_label = label_cells(data_run,cluster_IDs)
    gene_list = c('GZMA','GZMB','GZMH','GZMK','FCER1G' ,'CXCR4','KLRF1', 
                  'KLRB1', 'KLRD1', 'KLRC1', 'KLRG1', 'IL2RB', 'IL2RG', 'TSC22D3', 'NR4A2', 
                  'EVL', 'IFITM2', 'TNFAIP3','TGFB1','NFKBIA', 'GNLY', 'NKG7','FCGR3A', 'CCND3'
                  , 'LTB','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
                  'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1')
    
    gene_list = c('KLRK1','ITGAX','CX3CR1','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
                  'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1')
    folder_name = 'Dexa Gene Feature Plots'

    FeaturePlot_GeneList(data_run_label,folder_base_output,folder_name)
      
      
    
    
    
    # 
    # gene_list = c('GZMA','GZMB','GZMH','GZMK','FCER1G' ,'CXCR4','KLRF1', 
    #               'KLRB1', 'KLRD1', 'KLRC1', 'KLRG1', 'IL2RB', 'IL2RG', 'TSC22D3', 'NR4A2', 
    #               'EVL', 'IFITM2', 'TNFAIP3','TGFB1','NFKBIA', 'GNLY', 'NKG7','FCGR3A', 'CCND3',
    #               'LTB','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4','GIMAP7',
    #               'KLRK1','ITGAX','CX3CR1','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
    #               'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1','IL2RB')
    # 
###############################################################################
    
    gene_list = c('CXCR4','NFKBIA', 'TNFAIP3', 'NR4A2', 'TGFB1', 'RGS1','TSC22D3') 
    
    # HeatMap
    for (category in category_list){
      folder_name = paste0('Dexa Gene DoHeatmap ',category)
      folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
      dir.create( folder_featureplot, recursive = TRUE)
      data = SubsetData(object = data_run, cells = data_run$split_var == category )
      print(category)
      print(data)
      data = subset(data, idents = c('CD8+ T Cell','NK','CD14+ Mono'))
      
      plot  = DoHeatmap(object = data, features = gene_list,
                        group.by = "ident", disp.min = 0, disp.max = 2.5) +
        ggtitle(category ) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(size=24))
      
      pathName <- paste0(folder_featureplot,paste0(category,' All','.png'))  
      png(file=pathName,width=600, height=600)
      print(plot)
      dev.off()
      
      # for (i in 1:(length(gene_list))){
      #   feature = as.character(gene_list[i])
      #   print(feature)
      #   
      #   plot  = DoHeatmap(object = data, features = feature,
      #           group.by = "ident") +
      #           ggtitle(category ) + 
      #           theme(plot.title = element_text(hjust = 0.5)) +
      #           theme(plot.title = element_text(size=24))
      #     
      #   pathName <- paste0(folder_featureplot,paste0( feature,' ',category,'.png'))  
      #   png(file=pathName,width=1500, height=600)
      #   print(plot)
      #   dev.off()
      # }
    }  
 #################################################   
    # HeatMap
    gene_list = c('HLA-DQB1', 'HLA-DRB5', 'FCGR1A', 'CCR2', 'CX3CR1') 
    for (category in category_list){
      folder_name = paste0('Dexa Gene DoHeatmap ',category)
      folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
      dir.create( folder_featureplot, recursive = TRUE)
      data = SubsetData(object = data_run, cells = data_run$split_var == category )
      data = subset(data, idents = c('CD14+ Mono'))
      
      plot  = DoHeatmap(object = data, features = gene_list,
                        group.by = "ident", disp.min = 0, disp.max = 2.5) +
        ggtitle(category ) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(size=24))
      
      pathName <- paste0(folder_featureplot,paste0(category,' All Mono','.png'))  
      png(file=pathName,width=600, height=600)
      print(plot)
      dev.off()
      
      for (i in 1:(length(gene_list))){
        feature = as.character(gene_list[i])
        print(feature)
        
        # plot  = DoHeatmap(object = data, features = feature,
        #                   group.by = "ident") +
        #   ggtitle(category ) + 
        #   theme(plot.title = element_text(hjust = 0.5)) +
        #   theme(plot.title = element_text(size=24))
        # 
        # pathName <- paste0(folder_featureplot,paste0( feature,' ',category,'.png'))  
        # png(file=pathName,width=1500, height=600)
        # print(plot)
        # dev.off()
      }
    }    
    
    
    ###############################
    ## Save Data in matrix
    ###############################
    
    #SaveAsMatrix(data_run,folder_base_output)
    
    ###############################
    ## Plot samples seperately
    ###############################
    sample_list = unique(data_run$sample_name)
    data_run_label = label_cells(data_run,cluster_IDs)
    
    for (sample in sample_list){
      folder_output = paste0(folder_base_output,'Samples Seperate/',sample,'/')
      print(folder_output)
      pathName <- paste0(folder_output,'Cluster')
      dir.create( pathName, recursive = TRUE)
    
      pathName <- paste0(folder_output,'DE')
      dir.create( pathName, recursive = TRUE)
      
      pathName <- paste0(folder_output,'Stats')
      dir.create( pathName, recursive = TRUE)
      
      pathName <- paste0(folder_output,'PCA')
      dir.create( pathName, recursive = TRUE)
      
      data = SubsetData(object = data_run_label, cells = data_run$sample_name == sample )
      
      plotAll(data,folder_output,sample_name,sampleParam,label_TF = TRUE,integrate_TF = FALSE,  DE_perm_TF = FALSE,clusterTF = FALSE)
    }
    
  }
}
