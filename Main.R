# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Main.R')

# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

run = FALSE

folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all


for(i in patient_list){
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  if (run == TRUE){
    print('Running')
    data_orig = load_data(filename)
    param_data = sampleParam[ sampleParam$Sample ==sample_name , ] 
    
    
    filter <- TRUE
    regress_TF = FALSE
    folder = makeFolders(folder_base,sample_name,filter,regress_TF,TRUE)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    
    save(data,file=paste0(folder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    
    get_cellType(data,data_orig,folder,sample_name)
  
    
    regress_TF = TRUE
    folder = makeFolders(folder_base,sample_name,filter,regress_TF,TRUE)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)

    save(data,file=paste0(folder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)

    get_cellType(data,data_orig,folder,sample_name)
    
    filter <- FALSE
    regress_TF = FALSE
    folder = makeFolders(folder_base,sample_name,filter,regress_TF,TRUE)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    
    save(data,file=paste0(folder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    
    get_cellType(data,data_orig,folder,sample_name)
    
    
    regress_TF = TRUE
    folder = makeFolders(folder_base,sample_name,filter,regress_TF,TRUE)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    
    save(data,file=paste0(folder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,'parameters.csv'),row.names=FALSE)
    
    get_cellType(data,data_orig,folder,sample_name)
    
    

    print('done!')
    
    
    
  }else{
    runAll = FALSE
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name))
    
    if (!runAll){
      filter = sampleParam$filter_TF[sampleParam['Sample'] == sample_name]
      regress_TF = sampleParam$regress_TF[sampleParam['Sample'] == sample_name]
    

      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
    }else{
      # Filter on
      filter = TRUE
      regress_TF = FALSE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE)
      print(paste0('folder: ', folder))
       
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      
      # Filter off
      filter = FALSE
      regress_TF = FALSE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      
      regress_TF = TRUE
      folder = makeFolders(folder_base,sample_name,filter,regress_TF,FALSE)
      print(paste0('folder: ', folder))
      
      data <- loadRData(paste0(folder,'data.Robj'))
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = FALSE)
      plotAll(data,folder,sample_name,filter,regress_TF, label_TF = TRUE)
      }
      #get_cellType(data,data_orig,folder,sample_name,filter)
  }
}

