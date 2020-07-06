# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Individual/Main_individual.R')
gc()

# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(stringr)

#library(biomaRt)

base = '/home/sujwary/Desktop/scRNA/'
source(paste0(base,'Code/Functions.R' ))
source(paste0(base,'Code/Plot_func.R'))
source(paste0(base,'Code/Integration/PlotAll.R'))
source(paste0(base,'Code/Individual/run_pipeline.R'))

run = F

folder_base_output <- '/home/sujwary/Desktop/scRNA/Output/'

filename_sampleParam <- paste0(base,'Data/sample_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- paste0(base,'Data/EloRD Meta.xlsx')
metaData <- read_excel(filename_metaData)
#metaData = metaData[metaData$`CHIP sample` != 'EOT',]

cell_features = getCellMarkers(base)
i = 62
for(i in 19){ #30:33
  sample_name <- metaData$Sample[i]
  filename <- paste(base,"Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  param_data = sampleParam[ sampleParam$Sample ==sample_name , ] 
  if (run == TRUE){
    #browser()
    print('Running')
    data_orig = load_data(filename)

    filter <- F
    regress_TF = F
    folder = makeFolders(folder_base_output,sample_name,filter,regress_TF, makeFolder_TF= TRUE,nFeature_RNA_list,percent_mt)
    filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
    data = run_pipeline(data_orig,folder_base_output,sample_name,sampleParam,filter,regress_TF)

    save(data,file=paste0(folder,paste0('data_run_PCAdim',PCA_dim,'.Robj')))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)

    PlotKnownMarkers(data, paste0(filepath_cluster,'Cell Type/Heatmap/'), cell_features = cell_features,
                     plotType ='HeatMap' ,split = FALSE,featurePlotFix = TRUE, str = '')
    PlotKnownMarkers(data, paste0(filepath_cluster,'Cell Type/FeaturePlot/'), cell_features = cell_features,
                     plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = FALSE, str = '')
    
    filter = T
    regress_TF = F
    folder = makeFolders(folder_base_output,sample_name,filter,regress_TF, makeFolder_TF= TRUE,nFeature_RNA_list,percent_mt)
    filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
    data = run_pipeline(data_orig,folder_base_output,sample_name,sampleParam,filter,regress_TF)

    save(data,file=paste0(folder,paste0('data_run_PCAdim',PCA_dim,'.Robj')))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)

    PlotKnownMarkers(data, paste0(filepath_cluster,'Cell Type/Heatmap/'), cell_features = cell_features,
                     plotType ='HeatMap' ,split = FALSE,featurePlotFix = TRUE, str = '')
    PlotKnownMarkers(data, paste0(filepath_cluster,'Cell Type/FeaturePlot/'), cell_features = cell_features,
                     plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = FALSE, str = '')


  

    print('done!')
    
    
    
  }else{
    runAll = FALSE
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name))
    
    if (!runAll){
      #filter = sampleParam$filter_TF[sampleParam['Sample'] == sample_name]
      #regress_TF = sampleParam$regress_TF[sampleParam['Sample'] == sample_name]
      resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
      PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
      
      folder = makeFolders(folder_base_output,sample_name,filter = T,regress_TF = T, makeFolder_TF= TRUE,nFeature_RNA_list,percent_mt)      
      print(paste0('folder: ', folder))
      
      data = loadRData(file=paste0(folder,paste0('data_run_PCAdim',PCA_dim,'.Robj')))
      data$split_var = ''
      data$Response = ''
      plotAll(data,folder,sample_name,sampleParam,label_TF = F, cell_features = cell_features,integrate_TF = FALSE, DE_perm_TF = FALSE,
              clusterTF = T, markersTF = TRUE,keepOldLabels = F,groupBy = NA,
              PCA_dim = NA, resolution_val = NA)
      
      
      plotAll(data,folder,sample_name,sampleParam,label_TF = T, cell_features = cell_features,integrate_TF = FALSE, DE_perm_TF = FALSE,
              clusterTF = T, markersTF = TRUE,keepOldLabels = F,groupBy = NA,
              PCA_dim = NA, resolution_val = NA)
      

      data = getCluster (data,resolution_val, PCA_dim)
      filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
      write.csv(param_data, file = paste0(filepath_cluster,'parameters.csv'),row.names=FALSE)
      
      PlotKnownMarkers(data,data, paste0(filepath_cluster,'Cell Type/Heatmap/'), cell_features = NA,
                       plotType ='HeatMap' ,split = FALSE,featurePlotFix = TRUE, str = '')
      PlotKnownMarkers(data,data, paste0(filepath_cluster,'Cell Type/FeaturePlot/'), cell_features = NA,
                       plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = FALSE, str = '')

    }else{
      # Filter on
      filter = TRUE
      regress_TF = FALSE
      folder = makeFolders(folder_base_output,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
       
      data = loadRData(paste0(folder,paste0('data_run_PCAdim',PCA_dim,'.Robj')))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base_output,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      
      # Filter off
      filter = FALSE
      regress_TF = FALSE
      folder = makeFolders(folder_base_output,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base_output,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      }
      #get_cellType(data,data_orig,folder,cell_features = NA)
  }
}

