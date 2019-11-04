
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

run = TRUE

folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

for(i in 20){
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  if (run == TRUE){
    print('running')
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

    print('done!')
    
    
    
  }else{
    
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name))
    
    
    filter = sampleParam$filter_TF[sampleParam['Sample'] == sample_name]
    regress_TF = sampleParam$regress_TF[sampleParam['Sample'] == sample_name]
    
    folder = makeFolders(folder_base,sample_name,filter,regress_TF,TRUE)
    print(paste0('folder: ', folder))
     
    data <- loadRData(paste0(folder,'data.Robj'))
    label_TF = FALSE
    plotAll(data,folder,sample_name,filter,regress_TF, label_TF)
    label_TF = TRUE
    plotAll(data,folder,sample_name,filter,regress_TF, label_TF)
    #get_cellType(data,data_orig,folder,sample_name,filter)
  }
}

