# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Main.R')
gc()

# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Individual/run_pipeline.R')

run = TRUE

folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

for(i in 4){ #30:33
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  if (run == TRUE){
    #browser()
    print('Running')
    data_orig = load_data(filename)
    param_data = sampleParam[ sampleParam$Sample ==sample_name , ] 
    
    filter <- FALSE
    regress_TF = FALSE
    folder = makeFolders(folder_base,sample_name,filter,regress_TF, makeFolder_TF = TRUE)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
     
    save(data,file=paste0(folder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
     
    MakeFeaturePlot(data,data_orig,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)    
    
    
    # filter <- TRUE
    # regress_TF = FALSE
    # folder = makeFolders(folder_base,sample_name,filter,regress_TF, makeFolder_TF= TRUE)
    # data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    # 
    # 
    # save(data,file=paste0(folder,'data.Robj'))
    # write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    # 
    # MakeFeaturePlot(data,data_orig,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)
    # 
    # regress_TF = TRUE
    # folder = makeFolders(folder_base,sample_name,filter,regress_TF, makeFolder_TF = TRUE)
    # data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    # 
    # save(data,file=paste0(folder,'data.Robj'))
    # write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    # 
    # MakeFeaturePlot(data,data_orig,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)
    # 
    # filter <- FALSE
    # regress_TF = FALSE
    # folder = makeFolders(folder_base,sample_name,filter,regress_TF, makeFolder_TF = TRUE)
    # data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    # 
    # save(data,file=paste0(folder,'data.Robj'))
    # write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    # 
    # MakeFeaturePlot(data,data_orig,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)    
    # 
    # regress_TF = TRUE
    # folder = makeFolders(folder_base,sample_name,filter,regress_TF, makeFolder_TF = TRUE)
    # data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    # 
    # save(data,file=paste0(folder,'data.Robj'))
    # write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    # 
    # MakeFeaturePlot(data,data_orig,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)    
    

    print('done!')
    
    
    
  }else{
    runAll = FALSE
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name))
    
    if (!runAll){
      filter = sampleParam$filter_TF[sampleParam['Sample'] == sample_name]
      regress_TF = sampleParam$regress_TF[sampleParam['Sample'] == sample_name]
      resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
      PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
      
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data = loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,sampleParam,label_TF = FALSE,integrate_TF = FALSE, DE_perm_TF = TRUE)
      plotAll(data,folder,sample_name,sampleParam,label_TF = TRUE,integrate_TF = FALSE, DE_perm_TF = TRUE)
      
      
      data = getCluster (data,resolution_val, PCA_dim)
      MakeFeaturePlot(data,data,paste0(folder,'Cell Type/'), cell_features = NA, split = FALSE)

    }else{
      # Filter on
      filter = TRUE
      regress_TF = FALSE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
       
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      
      # Filter off
      filter = FALSE
      regress_TF = FALSE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE, makeFolder_TF = TRUE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      }
      #get_cellType(data,data_orig,folder,cell_features = NA)
  }
}

